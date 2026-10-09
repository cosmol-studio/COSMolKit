// Reference preparation only; never built or invoked by ordinary Cargo tests.
// Compile with the pinned RDKit SDK/libraries and pass torsion.etkdg.v2.mol.
// Calls native algorithms directly; no CK outputs or implementation are used.
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/versions.h>
#include <iomanip>
#include <iostream>
#include <memory>

void points(const char *key, const RDGeom::PointPtrVect &positions) {
  std::cout << "\"" << key << "\":[";
  for (unsigned i = 0; i < positions.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << '[' << (*positions[i])[0] << ',' << (*positions[i])[1]
              << ',' << (*positions[i])[2] << ']';
  }
  std::cout << "],";
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  if (std::string(RDKit::rdkitVersion) != "2026.03.6") return 4;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::MolFileToMol(argv[1], true, false));
  if (!mol) return 3;
  const auto count = mol->getNumAtoms();
  DistGeom::BoundsMatPtr bounds(new DistGeom::BoundsMatrix(count));
  RDKit::DGeomHelpers::initBoundsMat(bounds);
  RDKit::DGeomHelpers::setTopolBounds(*mol, bounds, true, false, false, true);
  const bool smoothing = DistGeom::triangleSmoothBounds(bounds);
  std::vector<RDGeom::Point3D> storage(count);
  RDGeom::PointPtrVect positions;
  RDGeom::Point3DPtrVect positions3d;
  for (auto &point : storage) {
    positions.push_back(&point);
    positions3d.push_back(&point);
  }
  RDKit::rng_type generator(42);
  RDKit::uniform_double distribution(0., 1.);
  RDKit::double_source_type rng(generator, distribution);
  RDNumeric::SymmMatrix<double> distances(count);
  DistGeom::pickRandomDistMat(*bounds, distances, rng);
  bool initial_ok = DistGeom::computeInitialCoords(distances, positions, rng, true, 1);
  std::cout << std::setprecision(17) << "{\"smoothing\":" << smoothing
            << ",\"initial_ok\":" << initial_ok << ',';
  points("initial", positions);
  boost::dynamic_bitset<> fixed(count);
  std::unique_ptr<ForceFields::ForceField> first(DistGeom::constructForceField(
      *bounds, positions, DistGeom::VECT_CHIRALSET{}, 1., .1, nullptr, 5., &fixed));
  first->initialize();
  double e1 = first->calcEnergy();
  if (e1 > 1e-5) while (first->minimize(400, 1e-3, 1e-6) != 0) {}
  double e2 = first->calcEnergy();
  points("first_minimized", positions);
  ForceFields::CrystalFF::CrystalFFDetails details{};
  details.boundsMatForceScaling = 1.;
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, false, false, 2, false);
  for (const auto atom : mol->atoms()) details.atomNums.push_back(atom->getAtomicNum());
  DistGeom::BoundsMatPtr topology_bounds(new DistGeom::BoundsMatrix(count));
  RDKit::DGeomHelpers::initBoundsMat(topology_bounds);
  RDKit::DGeomHelpers::setTopolBounds(*mol, topology_bounds, details.bonds, details.angles,
                                     true, false, false, true, true, true);
  std::unique_ptr<ForceFields::ForceField> second(DistGeom::construct3DForceField(*bounds, positions3d, details));
  second->initialize();
  double e3 = second->calcEnergy();
  if (e3 > 1e-5) second->minimize(300, 1e-3, 1e-6);
  double e4 = second->calcEnergy();
  points("final", positions);
  std::cout << "\"first_energy_before\":" << e1 << ",\"first_energy_after\":" << e2
            << ",\"torsion_energy_before\":" << e3 << ",\"torsion_energy_after\":" << e4 << "}\n";
}
