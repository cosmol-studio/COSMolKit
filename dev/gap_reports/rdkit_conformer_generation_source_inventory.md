# RDKit Conformer Generation Source Inventory

This inventory is the Step 3 artifact for
`dev/rdkit_conformer_generation_full_port_plan.md`.

The implementation target is full source-backed RDKit conformer generation
coverage for `EmbedMolecule` and `EmbedMultipleConfs`. This is an inventory
only; reuse eligibility is handled separately by
`rdkit_conformer_generation_reuse_map.md`.

## Baseline Source Files

The conformer-generation baseline includes the following RDKit files:

- `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`
- `third_party/rdkit/Code/DistGeom/ChiralSet.h`
- `third_party/rdkit/Code/DistGeom/ChiralViolationContribs.h`
- `third_party/rdkit/Code/DistGeom/ChiralViolationContribs.cpp`
- `third_party/rdkit/Code/DistGeom/DistViolationContribs.h`
- `third_party/rdkit/Code/DistGeom/DistViolationContribs.cpp`
- `third_party/rdkit/Code/DistGeom/FourthDimContribs.h`
- `third_party/rdkit/Code/DistGeom/DistGeomUtils.h`
- `third_party/rdkit/Code/DistGeom/DistGeomUtils.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Embedder.h`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Embedder.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/EmbedderUtils.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_v1.in`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_v2.in`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_smallrings.in`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_macrocycles.in`

## DistGeom Matrix And Smoothing Symbols

These symbols are already part of the conformer-generation baseline because
`EmbedMolecule` builds, smooths, and consumes bounds matrices.

| Source | Required symbol |
|---|---|
| `BoundsMatrix.h` | `DistGeom::BoundsMatrix` |
| `BoundsMatrix.h` | `DistGeom::BoundsMatPtr` |
| `BoundsMatrix.h` | `BoundsMatrix::BoundsMatrix(unsigned int)` |
| `BoundsMatrix.h` | `BoundsMatrix::getData()` |
| `BoundsMatrix.h` | `BoundsMatrix::numRows()` |
| `BoundsMatrix.h` | `BoundsMatrix::getUpperBound(unsigned int, unsigned int)` |
| `BoundsMatrix.h` | `BoundsMatrix::setUpperBound(unsigned int, unsigned int, double)` |
| `BoundsMatrix.h` | `BoundsMatrix::setUpperBoundIfBetter(unsigned int, unsigned int, double)` |
| `BoundsMatrix.h` | `BoundsMatrix::getLowerBound(unsigned int, unsigned int)` |
| `BoundsMatrix.h` | `BoundsMatrix::setLowerBound(unsigned int, unsigned int, double)` |
| `BoundsMatrix.h` | `BoundsMatrix::setLowerBoundIfBetter(unsigned int, unsigned int, double)` |
| `BoundsMatrix.h` | `BoundsMatrix::checkValid()` |
| `TriangleSmooth.h/.cpp` | `DistGeom::triangleSmoothBounds(BoundsMatrix*, double)` |
| `TriangleSmooth.h/.cpp` | `DistGeom::triangleSmoothBounds(BoundsMatPtr, double)` |

## DistGeom Chiral And Violation-Contribution Symbols

These symbols are required by the DG cleanup force field, tetrahedral checks,
chiral checks, and fourth-dimensional chirality minimization.

| Source | Required symbol |
|---|---|
| `ChiralSet.h` | `DistGeom::ChiralSetStructureFlags` |
| `ChiralSet.h` | `DistGeom::ChiralSet` |
| `ChiralSet.h` | `DistGeom::ChiralSetPtr` |
| `ChiralSet.h` | `DistGeom::VECT_CHIRALSET` |
| `ChiralSet.h` | `ChiralSet::ChiralSet(unsigned int, unsigned int, unsigned int, unsigned int, int, double, std::uint64_t)` |
| `ChiralViolationContribs.h/.cpp` | `DistGeom::calcChiralVolume(unsigned int, unsigned int, unsigned int, unsigned int, double*)` |
| `ChiralViolationContribs.h/.cpp` | `DistGeom::calcChiralVolume(unsigned int, unsigned int, unsigned int, unsigned int, RDGeom::PointPtrVect&)` |
| `ChiralViolationContribs.h` | `DistGeom::ChiralViolationContribsParams` |
| `ChiralViolationContribs.h/.cpp` | `DistGeom::ChiralViolationContribs` |
| `ChiralViolationContribs.cpp` | `ChiralViolationContribs::ChiralViolationContribs(ForceField*)` |
| `ChiralViolationContribs.cpp` | `ChiralViolationContribs::addContrib(const ChiralSet*, double)` |
| `ChiralViolationContribs.cpp` | `ChiralViolationContribs::getEnergy(double*) const` |
| `ChiralViolationContribs.cpp` | `ChiralViolationContribs::getGrad(double*, double*) const` |
| `DistViolationContribs.h` | `DistGeom::DistViolationContribsParams` |
| `DistViolationContribs.h/.cpp` | `DistGeom::DistViolationContribs` |
| `DistViolationContribs.cpp` | local helper `distance2(unsigned int, unsigned int, double*)` |
| `DistViolationContribs.cpp` | local helper `distance(unsigned int, unsigned int, double*)` |
| `DistViolationContribs.h` | `DistViolationContribs::addContrib(unsigned int, unsigned int, double, double, double)` |
| `DistViolationContribs.cpp` | `DistViolationContribs::getEnergy(double*) const` |
| `DistViolationContribs.cpp` | `DistViolationContribs::getGrad(double*, double*) const` |
| `FourthDimContribs.h` | `DistGeom::FourthDimContribsParams` |
| `FourthDimContribs.h` | `DistGeom::FourthDimContribs` |
| `FourthDimContribs.h` | `FourthDimContribs::addContrib(unsigned int, double)` |
| `FourthDimContribs.h` | `FourthDimContribs::getEnergy(double*) const` |
| `FourthDimContribs.h` | `FourthDimContribs::getGrad(double*, double*) const` |

## DistGeom Utility Symbols

These symbols implement random distance-matrix generation, metric embedding,
random box coordinates, and DG/ETKDG cleanup force fields.

| Source | Required symbol |
|---|---|
| `DistGeomUtils.cpp` | constant `EIGVAL_TOL` |
| `DistGeomUtils.cpp` | constant `KNOWN_DIST_TOL` |
| `DistGeomUtils.cpp` | constant `KNOWN_DIST_FORCE_CONSTANT` |
| `DistGeomUtils.h/.cpp` | `DistGeom::pickRandomDistMat(const BoundsMatrix&, SymmMatrix<double>&, int)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::pickRandomDistMat(const BoundsMatrix&, SymmMatrix<double>&, RDKit::double_source_type&)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::computeInitialCoords(const SymmMatrix<double>&, PointPtrVect&, bool, unsigned int, int)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::computeInitialCoords(const SymmMatrix<double>&, PointPtrVect&, RDKit::double_source_type&, bool, unsigned int)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::computeRandomCoords(PointPtrVect&, double, int)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::computeRandomCoords(PointPtrVect&, double, RDKit::double_source_type&)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::constructForceField(const BoundsMatrix&, PointPtrVect&, const VECT_CHIRALSET&, double, double, map<pair<int,int>,double>*, double, dynamic_bitset<>*)` |
| `DistGeomUtils.cpp` | `DistGeom::addImproperTorsionTerms(ForceField*, const vector<vector<int>>&, dynamic_bitset<>&)` |
| `DistGeomUtils.cpp` | `DistGeom::addExperimentalTorsionTerms(ForceField*, const CrystalFFDetails&, const BoundsMatrix&, Point3DPtrVect&, double)` |
| `DistGeomUtils.cpp` | `DistGeom::add12Terms(ForceField*, const CrystalFFDetails&, const BoundsMatrix&, Point3DPtrVect&, dynamic_bitset<>&, double)` |
| `DistGeomUtils.cpp` | `DistGeom::add13Terms(ForceField*, const CrystalFFDetails&, const dynamic_bitset<>&, const dynamic_bitset<>&, bool, const BoundsMatrix&, Point3DPtrVect&, dynamic_bitset<>&, double)` |
| `DistGeomUtils.cpp` | `DistGeom::addLongRangeDistanceConstraints(ForceField*, const CrystalFFDetails&, const dynamic_bitset<>&, Point3DPtrVect&, double, const BoundsMatrix&, double)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::construct3DForceField(const BoundsMatrix&, Point3DPtrVect&, const CrystalFFDetails&)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::construct3DForceField(const BoundsMatrix&, Point3DPtrVect&, const CrystalFFDetails&, const map<pair<unsigned int,unsigned int>,double>&)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::constructPlain3DForceField(const BoundsMatrix&, Point3DPtrVect&, const CrystalFFDetails&)` |
| `DistGeomUtils.h/.cpp` | `DistGeom::construct3DImproperForceField(const BoundsMatrix&, Point3DPtrVect&, const vector<vector<int>>&, const vector<vector<int>>&, const vector<int>&)` |
| `DistGeomUtils.h` | inline `DistGeom::construct3DImproperForceField(const BoundsMatrix&, Point3DPtrVect&, const CrystalFFDetails&)` |

## BoundsMatrixBuilder Symbols

These symbols are required by `setupInitialBoundsMatrix` and by wrapper
`GetMoleculeBoundsMatrix`. Some are already within the existing DG bounds
port, but they remain part of this conformer-generation baseline until Step 5
proves reuse.

| Source | Required symbol |
|---|---|
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::Path14Configuration` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::ComputedData` |
| `BoundsMatrixBuilder.cpp` | `ComputedData::visitedBound` |
| `BoundsMatrixBuilder.cpp` | `ComputedData::setTopolBounds` state fields |
| `BoundsMatrixBuilder.h/.cpp` | `DGeomHelpers::initBoundsMat(BoundsMatrix*, double, double)` |
| `BoundsMatrixBuilder.h/.cpp` | `DGeomHelpers::initBoundsMat(BoundsMatPtr, double, double)` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::set12Bounds` overloads |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::set13Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::set14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::set15Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::setLowerBoundVDW` overloads |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkAndSetBounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_set13BoundsHelper` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setRingAngle` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_getAtomStereo` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setInRing14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setTwoInSameRing14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setTwoInDiffRing14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setShareRingBond14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkH2NX3H1OX2` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkNhChChNh` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkAmideEster14` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_isCarbonyl` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkAmideEster15` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setChain14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_record14Path` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_setMacrocycleAllInSameRing14Bounds` |
| `BoundsMatrixBuilder.h/.cpp` | `DGeomHelpers::setTopolBounds` overloads |
| `BoundsMatrixBuilder.h/.cpp` | `DGeomHelpers::collectBondsAndAngles` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_compute15DistsCisCis` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_compute15DistsCisTrans` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_compute15DistsTransTrans` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_compute15DistsTransCis` |
| `BoundsMatrixBuilder.cpp` | `DGeomHelpers::_set15BoundsHelper` |

## Embedder Types, Constants, And Presets

These symbols define the public and internal controls for RDKit embedding.

| Source | Required symbol |
|---|---|
| `Embedder.h` | `DGeomHelpers::EmbedFailureCauses` |
| `Embedder.h` | `DGeomHelpers::EmbedParameters` |
| `Embedder.h` | `EmbedParameters::maxIterations` |
| `Embedder.h` | `EmbedParameters::numThreads` |
| `Embedder.h` | `EmbedParameters::randomSeed` |
| `Embedder.h` | `EmbedParameters::clearConfs` |
| `Embedder.h` | `EmbedParameters::useRandomCoords` |
| `Embedder.h` | `EmbedParameters::boxSizeMult` |
| `Embedder.h` | `EmbedParameters::randNegEig` |
| `Embedder.h` | `EmbedParameters::numZeroFail` |
| `Embedder.h` | `EmbedParameters::coordMap` |
| `Embedder.h` | `EmbedParameters::optimizerForceTol` |
| `Embedder.h` | `EmbedParameters::ignoreSmoothingFailures` |
| `Embedder.h` | `EmbedParameters::enforceChirality` |
| `Embedder.h` | `EmbedParameters::useExpTorsionAnglePrefs` |
| `Embedder.h` | `EmbedParameters::useBasicKnowledge` |
| `Embedder.h` | `EmbedParameters::verbose` |
| `Embedder.h` | `EmbedParameters::basinThresh` |
| `Embedder.h` | `EmbedParameters::pruneRmsThresh` |
| `Embedder.h` | `EmbedParameters::onlyHeavyAtomsForRMS` |
| `Embedder.h` | `EmbedParameters::ETversion` |
| `Embedder.h` | `EmbedParameters::boundsMat` |
| `Embedder.h` | `EmbedParameters::embedFragmentsSeparately` |
| `Embedder.h` | `EmbedParameters::useSmallRingTorsions` |
| `Embedder.h` | `EmbedParameters::useMacrocycleTorsions` |
| `Embedder.h` | `EmbedParameters::useMacrocycle14config` |
| `Embedder.h` | `EmbedParameters::timeout` |
| `Embedder.h` | `EmbedParameters::CPCI` |
| `Embedder.h` | `EmbedParameters::callback` |
| `Embedder.h` | `EmbedParameters::forceTransAmides` |
| `Embedder.h` | `EmbedParameters::useSymmetryForPruning` |
| `Embedder.h` | `EmbedParameters::boundsMatForceScaling` |
| `Embedder.h` | `EmbedParameters::trackFailures` |
| `Embedder.h` | `EmbedParameters::failures` |
| `Embedder.h` | `EmbedParameters::enableSequentialRandomSeeds` |
| `Embedder.h` | `EmbedParameters::symmetrizeConjugatedTerminalGroupsForPruning` |
| `Embedder.cpp` | constants `INTERRUPT_MESSAGE`, `M_PI_2`, `ERROR_TOL`, `MAX_MINIMIZED_E_PER_ATOM`, `MIN_TETRAHEDRAL_CHIRAL_VOL`, `TETRAHEDRAL_CENTERINVOLUME_TOL` |
| `Embedder.cpp` | preset `DGeomHelpers::KDG` |
| `Embedder.cpp` | preset `DGeomHelpers::ETDG` |
| `Embedder.cpp` | preset `DGeomHelpers::ETDGv2` |
| `Embedder.cpp` | preset `DGeomHelpers::ETKDG` |
| `Embedder.cpp` | preset `DGeomHelpers::ETKDGv2` |
| `Embedder.cpp` | preset `DGeomHelpers::ETKDGv3` |
| `Embedder.cpp` | preset `DGeomHelpers::srETKDGv3` |
| `Embedder.cpp` | `DGeomHelpers::detail::EmbedArgs` |

## Embedder Functions

These symbols are the main conformer-generation control flow and must be
source-reproduced without replacing behavior with heuristic approximations.

| Source | Required symbol |
|---|---|
| `Embedder.cpp` | local helper `haveOppositeSign(double, double)` |
| `Embedder.cpp` | local helper `failmutex_get()` |
| `Embedder.cpp` | local helper `failmutex_create()` |
| `Embedder.cpp` | local helper `GetFailMutex()` |
| `Embedder.cpp` | `EmbeddingOps::_volumeTest` |
| `Embedder.cpp` | `EmbeddingOps::_sameSide` |
| `Embedder.cpp` | `EmbeddingOps::_centerInVolume` atom-index overload |
| `Embedder.cpp` | `EmbeddingOps::_centerInVolume` chiral-set overload |
| `Embedder.cpp` | `EmbeddingOps::_boundsFulfilled` |
| `Embedder.cpp` | `EmbeddingOps::generateInitialCoords` |
| `Embedder.cpp` | `EmbeddingOps::firstMinimization` |
| `Embedder.cpp` | `EmbeddingOps::checkTetrahedralCenters` |
| `Embedder.cpp` | `EmbeddingOps::checkChiralCenters` |
| `Embedder.cpp` | `EmbeddingOps::minimizeFourthDimension` |
| `Embedder.cpp` | `EmbeddingOps::minimizeWithExpTorsions` |
| `Embedder.cpp` | `EmbeddingOps::doubleBondGeometryChecks` |
| `Embedder.cpp` | `EmbeddingOps::doubleBondStereoChecks` |
| `Embedder.cpp` | `EmbeddingOps::finalChiralChecks` |
| `Embedder.cpp` | `EmbeddingOps::embedPoints` |
| `Embedder.cpp` | `DGeomHelpers::findDoubleBonds` |
| `Embedder.cpp` | `DGeomHelpers::findChiralSets` |
| `Embedder.cpp` | `DGeomHelpers::adjustBoundsMatFromCoordMap` |
| `Embedder.cpp` | `DGeomHelpers::initETKDG` |
| `Embedder.cpp` | `DGeomHelpers::setupInitialBoundsMatrix` |
| `Embedder.cpp` | local helper `_fillAtomPositions` |
| `Embedder.cpp` | local helper `_isConfFarFromRest` |
| `Embedder.cpp` | template helper `multiplication_overflows_` |
| `Embedder.cpp` | local helper `embedHelper_` |
| `Embedder.cpp` | local helper `getMolSelfMatches` |
| `Embedder.cpp` | `DGeomHelpers::EmbedMultipleConfs(ROMol&, INT_VECT&, unsigned int, EmbedParameters&)` |
| `Embedder.h` | inline `DGeomHelpers::EmbedMultipleConfs(ROMol&, unsigned int, EmbedParameters&)` |
| `Embedder.h` | inline `DGeomHelpers::EmbedMolecule(ROMol&, EmbedParameters&)` |
| `Embedder.h` | legacy inline `DGeomHelpers::EmbedMolecule` overload with scalar parameters |
| `Embedder.h` | legacy inline `DGeomHelpers::EmbedMultipleConfs` overload writing an `INT_VECT&` |
| `Embedder.h` | legacy inline `DGeomHelpers::EmbedMultipleConfs` overload returning an `INT_VECT` |
| `EmbedderUtils.cpp` | `DGeomHelpers::updateEmbedParametersFromJSON` |
| `EmbedderUtils.cpp` | `DGeomHelpers::embedParametersToJSON` |

## rdDistGeom Wrapper Inventory

These wrappers define public Python semantics in RDKit and must be represented
in COSMolKit Python APIs where the underlying COSMolKit state model supports
the behavior.

| Source | Required wrapper symbol |
|---|---|
| `rdDistGeom.cpp` | `RDKit::PyEmbedParameters` |
| `rdDistGeom.cpp` | `PyEmbedParameters::getFailureCounts()` |
| `rdDistGeom.cpp` | `PyEmbedParameters::setBoundsMatrix()` |
| `rdDistGeom.cpp` | `PyEmbedParameters::setCPCI()` |
| `rdDistGeom.cpp` | `PyEmbedParameters::setCoordMap()` |
| `rdDistGeom.cpp` | `RDKit::EmbedMolecule` scalar-argument wrapper |
| `rdDistGeom.cpp` | `RDKit::EmbedMolecule2` params wrapper |
| `rdDistGeom.cpp` | `RDKit::EmbedMultipleConfs` scalar-argument wrapper |
| `rdDistGeom.cpp` | `RDKit::EmbedMultipleConfs2` params wrapper |
| `rdDistGeom.cpp` | `RDKit::getMolBoundsMatrix` |
| `rdDistGeom.cpp` | `RDKit::getETKDG` |
| `rdDistGeom.cpp` | `RDKit::getETKDGv2` |
| `rdDistGeom.cpp` | `RDKit::getETKDGv3` |
| `rdDistGeom.cpp` | `RDKit::getsrETKDGv3` |
| `rdDistGeom.cpp` | `RDKit::getKDG` |
| `rdDistGeom.cpp` | `RDKit::getETDG` |
| `rdDistGeom.cpp` | `RDKit::getETDGv2` |
| `rdDistGeom.cpp` | `RDKit::getExpTorsHelper` |
| `rdDistGeom.cpp` | `RDKit::getExpTorsHelperWithParams` |
| `rdDistGeom.cpp` | `RDKit::embedParametersToJSONHelper` |
| `rdDistGeom.cpp` | Boost.Python `EmbedParameters` class field exports |
| `rdDistGeom.cpp` | Boost.Python `EmbedMolecule`, `EmbedMultipleConfs`, `GetMoleculeBoundsMatrix`, parameter-factory, torsion-helper, and JSON exports |

## CrystalFF And ETKDG Dependency Inventory

These symbols are reachable through `initETKDG`, `minimizeWithExpTorsions`,
and `DistGeomUtils` construction of ETKDG force fields.

| Source | Required symbol |
|---|---|
| `TorsionPreferences.h` | `ForceFields::CrystalFF::ExpTorsionAngle` |
| `TorsionPreferences.h` | `ForceFields::CrystalFF::CrystalFFDetails` |
| `TorsionPreferences.cpp` | `ForceFields::CrystalFF::ExpTorsionAngleCollection` |
| `TorsionPreferences.cpp` | `ExpTorsionAngleCollection::getParams` |
| `TorsionPreferences.cpp` | `ExpTorsionAngleCollection::ExpTorsionAngleCollection` |
| `TorsionPreferences.h/.cpp` | `ForceFields::CrystalFF::getExperimentalTorsions(const ROMol&, CrystalFFDetails&, vector<pair<...>>&, bool, bool, bool, bool, unsigned int, bool)` |
| `TorsionPreferences.h/.cpp` | `ForceFields::CrystalFF::getExperimentalTorsions(const ROMol&, CrystalFFDetails&, bool, bool, bool, bool, unsigned int, bool)` |
| `TorsionAngleContribs.h` | `ForceFields::CrystalFF::TorsionAngleContribsParams` |
| `TorsionAngleContribs.h/.cpp` | `ForceFields::CrystalFF::TorsionAngleContribs` |
| `TorsionAngleContribs.cpp` | `TorsionAngleContribs::TorsionAngleContribs(ForceField*)` |
| `TorsionAngleContribs.cpp` | `TorsionAngleContribs::addContrib(unsigned int, unsigned int, unsigned int, unsigned int, const vector<double>&, const vector<int>&)` |
| `TorsionAngleContribs.cpp` | `TorsionAngleContribs::getEnergy(double*) const` |
| `TorsionAngleContribs.cpp` | `TorsionAngleContribs::getGrad(double*, double*) const` |
| `TorsionAngleContribs.h/.cpp` | `ForceFields::CrystalFF::calcTorsionEnergy(const vector<double>&, const vector<int>&, double)` |
| `TorsionAngleContribs.cpp` | local `calcTorsionEnergyM6(const vector<double>&, const vector<int>&, double)` |
| `TorsionAngleM6.h/.cpp` | `ForceFields::CrystalFF::TorsionAngleContribM6` |
| `TorsionAngleM6.cpp` | `TorsionAngleContribM6::TorsionAngleContribM6(ForceField*, unsigned int, unsigned int, unsigned int, unsigned int, const vector<double>&, const vector<int>&)` |
| `TorsionAngleM6.cpp` | `TorsionAngleContribM6::getEnergy(double*) const` |
| `TorsionAngleM6.cpp` | `TorsionAngleContribM6::getGrad(double*, double*) const` |
| `TorsionAngleM6.h/.cpp` | `ForceFields::CrystalFF::calcTorsionEnergyM6(const vector<double>&, const vector<int>&, double)` |
| `torsionPreferences_v1.in` | data table `torsionPreferencesV1` |
| `torsionPreferences_v2.in` | data table `torsionPreferencesV2` |
| `torsionPreferences_smallrings.in` | data table `torsionPreferencesSmallRings` |
| `torsionPreferences_macrocycles.in` | data table `torsionPreferencesMacrocycles` |

## Reachable Force-Field Core Dependencies

The ETKDG and DG cleanup force fields use RDKit force-field core and contrib
types. COSMolKit already has force-field modules, but Step 5 must prove exact
reuse before any of these are treated as covered.

- `ForceFields::ForceField`
- `ForceFields::ForceFieldContrib`
- `ForceFields::DistanceConstraintContrib`
- `ForceFields::AngleConstraintContrib`
- `ForceFields::UFF::InversionContrib`
- `ForceFields::CrystalFF::TorsionAngleContribs`
- `ForceFields::CrystalFF::TorsionAngleContribM6`
- `ForceFields::MMFF::Utils::calcTorsionGrad`
- `ForceFields::MMFF::Utils::calcTorsionCosPhi`
- `RDGeom::Point`, `RDGeom::Point3D`, `RDGeom::PointPtrVect`, `RDGeom::Point3DPtrVect`, and `RDGeom::Point3DConstPtrVect`
- `RDKit::double_source_type` and seeded random-source behavior
- `RDNumeric::SymmMatrix<double>`
- `RDNumeric::EigenSolvers::powerEigenSolver`
- `boost::dynamic_bitset<>`

## Required Test Sources

The RDKit tests and fixtures below are part of the source baseline for parity
and branch coverage.

| Source | Coverage role |
|---|---|
| `third_party/rdkit/Code/DistGeom/testDistGeom.cpp` | `pickRandomDistMat`, `computeInitialCoords`, low-level DG behavior |
| `third_party/rdkit/Code/GraphMol/DistGeomHelpers/testDgeomHelpers.cpp` | legacy scalar API, DG/ETDG/KDG/ETKDG variants, pruning, constrained embedding, UFF postchecks, many regression cases |
| `third_party/rdkit/Code/GraphMol/DistGeomHelpers/catch_tests.cpp` | modern ETKDGv3, macrocycle, stereo, symmetry, timeout, and issue-regression coverage |
| `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/testDistGeom.py` | Python wrapper behavior, `EmbedParameters`, JSON, factories, `GetMoleculeBoundsMatrix`, default ETKDG semantics |
| `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/testCrystalFF.cpp` | CrystalFF torsion-energy and torsion-preference coverage |

## Required Fixture Files

The following fixture files must be imported or reproduced with provenance for
COSMolKit parity tests:

- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.1.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.2.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.3.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3483968.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/atropisomers.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/chirality_failure_test.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/embedDistOpti.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/initCoords.random.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.new.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle1.etkdg.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol`

## Public Behavior Inventory

The full COSMolKit conformer-generation port must expose the following
behavior without RDKit interop:

- Single-conformer generation equivalent to RDKit `EmbedMolecule`.
- Multi-conformer generation equivalent to RDKit `EmbedMultipleConfs`.
- Parameter presets equivalent to RDKit `KDG`, `ETDG`, `ETDGv2`, `ETKDG`, `ETKDGv2`, `ETKDGv3`, and `srETKDGv3`.
- Full `EmbedParameters` behavior for seeds, attempts, random coordinates, eigenvalue handling, smoothing failures, chirality, ETKDG options, fragment handling, pruning, symmetry pruning, custom bounds matrices, coord maps, CPCI, callbacks where representable, timeout, failure tracking, and sequential random seeds.
- Value-style molecule transformation APIs that preserve the source molecule.
- Coordinate-domain operation metadata and invariant checks for generated conformer rows.
- Python APIs and stubs matching the completed Rust behavior.
- Native Rust and Python examples for conformer generation and force-field post-optimization.

## Initial Coverage Verdict

The inventory shows that COSMolKit currently has substantial DG bounds and
force-field code, including Rust files under:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/core.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/uff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/*`

However, this Step 3 inventory does not mark any of those reusable. Step 5
must compare each reused Rust symbol against the conformer-generation source
baseline before implementation can rely on it.
