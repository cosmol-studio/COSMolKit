#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

namespace RDKit::MMFF {
enum { CONSTANT = 1, DISTANCE = 2 };
}  // namespace RDKit::MMFF

namespace ForceFields::MMFF::Utils {
double calcVdWEnergy(const double dist, const double R_star_ij,
                     const double wellDepth) {
  double const vdw1 = 1.07;
  double const vdw1m1 = vdw1 - 1.0;
  double const vdw2 = 1.12;
  double const vdw2m1 = vdw2 - 1.0;
  double dist2 = dist * dist;
  double dist7 = dist2 * dist2 * dist2 * dist;
  double aTerm = vdw1 * R_star_ij / (dist + vdw1m1 * R_star_ij);
  double aTerm2 = aTerm * aTerm;
  double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
  double R_star_ij2 = R_star_ij * R_star_ij;
  double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
  double bTerm = vdw2 * R_star_ij7 / (dist7 + vdw2m1 * R_star_ij7) - 2.0;
  double res = wellDepth * aTerm7 * bTerm;

  return res;
}

double calcEleEnergy(unsigned int, unsigned int, double dist, double chargeTerm,
                     std::uint8_t dielModel, bool is1_4) {
  double corr_dist = dist + 0.05;
  double const diel = 332.0716;
  double const sc1_4 = 0.75;
  if (dielModel == RDKit::MMFF::DISTANCE) {
    corr_dist *= corr_dist;
  }
  return (diel * chargeTerm / corr_dist * (is1_4 ? sc1_4 : 1.0));
}
}  // namespace ForceFields::MMFF::Utils

static unsigned long long bits(double value) {
  return static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(value));
}

static void emitVdwRow(unsigned int row, double dist, double radius, double depth,
                       double output) {
  std::printf("VDW\t%02u\t%016llx\t%016llx\t%016llx\t%016llx\n", row,
              bits(dist), bits(radius), bits(depth), bits(output));
}

static void emitEleRow(unsigned int row, unsigned int idx1, unsigned int idx2,
                       double dist, double charge, std::uint8_t dielModel,
                       bool is1_4, double output) {
  std::printf("ELE\t%03u\t%08x\t%08x\t%016llx\t%016llx\t%02x\t%u\t%016llx\n",
              row, idx1, idx2, bits(dist), bits(charge),
              static_cast<unsigned int>(dielModel), is1_4 ? 1U : 0U,
              bits(output));
}

int main() {
  constexpr std::array<double, 5> vdwDistances{0.0, 0.125, 1.0, 3.5, 12.0};
  constexpr std::array<double, 2> vdwRadii{1.25, 4.0};
  constexpr std::array<double, 2> wellDepths{0.5, -0.0};
  unsigned int vdwCalls = 0;
  unsigned int vdwRows = 0;

  for (double dist : vdwDistances) {
    for (double radius : vdwRadii) {
      for (double depth : wellDepths) {
        for (unsigned int repeat = 0; repeat < 2; ++repeat) {
          const double output =
              ForceFields::MMFF::Utils::calcVdWEnergy(dist, radius, depth);
          ++vdwCalls;
          if (repeat == 0) {
            emitVdwRow(vdwRows++, dist, radius, depth, output);
          }
        }
      }
    }
  }

  constexpr double infinity = std::numeric_limits<double>::infinity();
  constexpr std::array<std::array<double, 3>, 4> vdwControls{{
      {{-0.0, 1.25, 0.5}},
      {{-0.0, 1.25, -0.0}},
      {{infinity, 1.25, 0.5}},
      {{infinity, 1.25, -0.0}},
  }};
  for (const auto &input : vdwControls) {
    const double output = ForceFields::MMFF::Utils::calcVdWEnergy(
        input[0], input[1], input[2]);
    ++vdwCalls;
    emitVdwRow(vdwRows++, input[0], input[1], input[2], output);
  }

  constexpr std::array<double, 3> eleDistances{0.0, 0.75, 8.0};
  constexpr std::array<double, 3> charges{-0.0, 0.625, -1.5};
  constexpr std::array<std::uint8_t, 4> dielModels{1, 2, 0, 255};
  constexpr std::array<bool, 2> oneFourFlags{false, true};
  constexpr unsigned int largestIndex =
      std::numeric_limits<unsigned int>::max();
  unsigned int eleCalls = 0;
  unsigned int eleRows = 0;

  for (double dist : eleDistances) {
    for (double charge : charges) {
      for (std::uint8_t dielModel : dielModels) {
        for (bool is1_4 : oneFourFlags) {
          for (unsigned int repeat = 0; repeat < 2; ++repeat) {
            const unsigned int idx1 = repeat == 0 ? 0U : largestIndex;
            const unsigned int idx2 = repeat == 0 ? 1U : largestIndex;
            const double output = ForceFields::MMFF::Utils::calcEleEnergy(
                idx1, idx2, dist, charge, dielModel, is1_4);
            ++eleCalls;
            if (repeat == 0) {
              emitEleRow(eleRows++, idx1, idx2, dist, charge, dielModel, is1_4,
                         output);
            }
          }
        }
      }
    }
  }

  for (std::uint8_t dielModel : dielModels) {
    for (bool is1_4 : oneFourFlags) {
      const double output = ForceFields::MMFF::Utils::calcEleEnergy(
          17U, 29U, 0.75, 0.0, dielModel, is1_4);
      ++eleCalls;
      emitEleRow(eleRows++, 17U, 29U, 0.75, 0.0, dielModel, is1_4, output);
    }
  }

  constexpr std::array<double, 2> boundaryDistances{-0.025, infinity};
  constexpr std::array<std::uint8_t, 2> boundaryModels{1, 2};
  for (double dist : boundaryDistances) {
    for (std::uint8_t dielModel : boundaryModels) {
      for (bool is1_4 : oneFourFlags) {
        const double output = ForceFields::MMFF::Utils::calcEleEnergy(
            17U, 29U, dist, 0.625, dielModel, is1_4);
        ++eleCalls;
        emitEleRow(eleRows++, 17U, 29U, dist, 0.625, dielModel, is1_4,
                   output);
      }
    }
  }

  std::printf("COUNTS\tvdw_calls=%u\tvdw_rows=%u\tele_calls=%u"
              "\tele_rows=%u\ttotal_calls=%u\n",
              vdwCalls, vdwRows, eleCalls, eleRows, vdwCalls + eleCalls);
  if (vdwCalls != 44 || vdwRows != 24 || eleCalls != 160 || eleRows != 88) {
    return 1;
  }
  return 0;
}
