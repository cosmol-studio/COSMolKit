// Test-only adapter. Numeric result leaves are unsigned binary64 bit patterns.
#include <GraphMol/Descriptors/GraphDescriptors.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#ifndef CK_EXACT_IPC
#include <GraphMol/Descriptors/BCUT.h>
#include <GraphMol/Descriptors/EState.h>
#include <GraphMol/Descriptors/MiscDescriptors.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#endif
#include <Eigen/Core>
#include <bit>
#include <boost/json/src.hpp>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace json = boost::json;
namespace desc = RDKit::Descriptors;

static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
static_assert(EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 4 &&
              EIGEN_MINOR_VERSION == 0);
#ifdef __FAST_MATH__
#error Native descriptor references require IEEE non-finite semantics.
#endif

namespace {
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

#ifdef CK_EXACT_IPC
bool matchesIpcEntropyProfile() {
  // K(2,11) and K(8,47) as disconnected matrix blocks have nonzero weights
  // 1, 398 and 8272. Resolve log at runtime: the same executable can load
  // a different libm implementation for this probability on another CPU.
  volatile double nearOne = 8272.0 / 8671.0;
  // The other branch uses table reduction. This coefficient probability of
  // disjoint stars of degrees 1, 2, 4, 4, 4, 4, 5, 5 distinguishes the
  // selected glibc result from a correctly rounded logarithm.
  volatile double tableReduced = 37120.0 / 135000.0;
  if (std::fegetround() != FE_TONEAREST ||
      bits(std::log(nearOne)) != 0xbfa81e83b1b71ba6ULL ||
      bits(std::log(tableReduced)) != 0xbff4a86c42c223a7ULL)
    return false;

  // Propane distinguishes separate multiply/add in the RDKit entropy sum
  // from FMA contraction in the translation unit that actually computes it.
  std::unique_ptr<RDKit::ROMol> propane(RDKit::SmilesToMol("CCC"));
  return propane &&
         bits(desc::calcIpc(*propane, true)) == 0x3fed62adf1ea257bULL;
}
#endif

json::array row(const std::vector<double> &values) {
  json::array result;
  for (double x : values)
    result.emplace_back(bits(x));
  return result;
}

json::value evaluate(const RDKit::ROMol &mol, const json::string &task,
                     const json::object &profile) {
  if (task == "ipc" || task == "avg_ipc") {
    std::vector<double> matrix;
    const auto &supplied = profile.at("distance_matrix");
    if (!supplied.is_null()) {
      const auto n = mol.getNumAtoms();
      const auto &rows = supplied.as_array();
      if (rows.size() != n)
        throw std::invalid_argument("distance matrix has " +
                                    std::to_string(rows.size()) +
                                    " rows; expected " + std::to_string(n));
      for (const auto &r : rows) {
        if (r.as_array().size() != n)
          throw std::invalid_argument("distance matrix column count");
        for (const auto &entry : r.as_array())
          matrix.push_back(entry.to_number<double>());
      }
    }
    return bits(desc::calcIpc(mol, task == "avg_ipc",
                              matrix.empty() ? nullptr : matrix.data()));
  }
#ifndef CK_EXACT_IPC
  if (task == "bcut2d")
    return row(desc::BCUT2D(mol));
  if (task == "estate_indices")
    return row(desc::calcEStateIndices(mol));
  if (task == "estate_vsa")
    return row(desc::calcEState_VSA(mol));
  if (task == "vsa_estate")
    return row(desc::calcVSA_EState(mol));
  if (task == "peoe_vsa")
    return row(desc::calcPEOE_VSA(mol));
  if (task == "estate_extrema") {
    const auto v = desc::calcEStates(mol);
    return row({v.minEState, v.maxEState, v.minAbsEState, v.maxAbsEState});
  }
  if (task == "partial_charge_extrema") {
    const auto v = desc::calcPartialCharges(mol);
    return row({v.minPartialCharge, v.maxPartialCharge, v.minAbsPartialCharge,
                v.maxAbsPartialCharge});
  }
  if (task == "gasteiger_charges") {
    RDKit::computeGasteigerCharges(mol,
                                   profile.at("iterations").to_number<int>(),
                                   profile.at("strict").as_bool());
    json::array result;
    for (const auto atom : mol.atoms()) {
      result.emplace_back(row(
          {atom->getProp<double>(RDKit::common_properties::_GasteigerCharge),
           atom->getProp<double>(
               RDKit::common_properties::_GasteigerHCharge)}));
    }
    return result;
  }
#endif
  throw std::invalid_argument("task incompatible with native reference");
}

json::object response(const json::object &request) {
  json::object result{{"id", request.at("id")}};
  std::unique_ptr<RDKit::RWMol> mol;
  try {
    mol.reset(
        RDKit::SmilesToMol(std::string(request.at("smiles").as_string())));
  } catch (const RDKit::MolSanitizeException &) {
    // Same input rejection as a null result from SmilesToMol.
  } catch (const RDKit::SmilesParseException &) {
  }
  if (!mol) {
    result["rejected"] = "parse_rejected";
    return result;
  }
  try {
    result["bits"] = evaluate(*mol, request.at("task").as_string(),
                              request.at("profile").as_object());
  } catch (const std::exception &error) {
    result["error"] = error.what();
  }
  return result;
}
} // namespace

int main(int argc, char **argv) {
#ifdef CK_EXACT_IPC
  if (!matchesIpcEntropyProfile()) {
    std::cerr << "exact-Ipc: incompatible entropy rounding (libm or FMA)\n";
    return 1;
  }
#endif
  if (argc == 2 && std::string_view(argv[1]) == "--build-info") {
    std::cout << json::object{{"protocol", 1},
                              {"kind", CK_REFERENCE_KIND},
                              {"revision", CK_RDKIT_REVISION},
                              {"adapter_sha256", CK_ADAPTER_SHA},
                              {"cmake_sha256", CK_CMAKE_SHA},
                              {"compiler", CK_CXX},
                              {"build_type", CK_BUILD},
                              {"flags", CK_CXX_FLAGS},
                              {"boost", BOOST_LIB_VERSION},
                              {"eigen", "3.4.0"},
                              {"simd", Eigen::SimdInstructionSetsInUse()},
                              {"rounding", std::fegetround()}}
              << '\n';
    return 0;
  }
  if (argc != 1)
    return 2;
  try {
    std::string line;
    while (std::getline(std::cin, line))
      std::cout << response(json::parse(line).as_object()) << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return std::cin.eof() && std::cout ? 0 : 1;
}
