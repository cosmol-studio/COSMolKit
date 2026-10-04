// Source-only literal driver for the pinned RDKit MMFF scalar energy owners.
// Pin: 351f8f378f8ad6bbd517980c38896e66bf907af8c
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ForceFields {
namespace MMFF {

constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180.0 / M_PI;
constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;

namespace Utils {

//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
// Source: Code/ForceField/MMFF/BondStretch.cpp:32-42
double calcBondStretchEnergy(const double r0, const double kb,
                             const double distance) {
  double distTerm = distance - r0;
  double distTerm2 = distTerm * distTerm;
  double const c1 = MDYNE_A_TO_KCAL_MOL;
  double const cs = -2.0;
  double const c3 = 7.0 / 12.0;

  return (0.5 * c1 * kb * distTerm2 *
          (1.0 + cs * distTerm + c3 * cs * cs * distTerm2));
}

//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
// Source: Code/ForceField/MMFF/AngleBend.cpp:43-57
double calcAngleBendEnergy(const double theta0, const double ka, bool isLinear,
                           const double cosTheta) {
  double angle = RAD2DEG * acos(cosTheta) - theta0;
  double const cb = -0.006981317;
  double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
  double res = 0.0;

  if (isLinear) {
    res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cosTheta);
  } else {
    res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
  }

  return res;
}

}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields

static std::uint64_t bits(double value) {
  return std::bit_cast<std::uint64_t>(value);
}

static void emitBond(unsigned int row, double r0, double kb, double distance,
                     double energy) {
  std::printf("B\t%03u\t%016llx\t%016llx\t%016llx\t%016llx\n", row,
              static_cast<unsigned long long>(bits(r0)),
              static_cast<unsigned long long>(bits(kb)),
              static_cast<unsigned long long>(bits(distance)),
              static_cast<unsigned long long>(bits(energy)));
}

static void emitAngle(unsigned int row, double theta0, double ka,
                      bool isLinear, double cosTheta, double energy) {
  std::printf("A\t%03u\t%016llx\t%016llx\t%u\t%016llx\t%016llx\n", row,
              static_cast<unsigned long long>(bits(theta0)),
              static_cast<unsigned long long>(bits(ka)), isLinear ? 1U : 0U,
              static_cast<unsigned long long>(bits(cosTheta)),
              static_cast<unsigned long long>(bits(energy)));
}

int main() {
  constexpr std::array<double, 3> r0Values{0.0, 1.25, 2.0};
  constexpr std::array<double, 3> kbValues{-0.0, 0.5, 4.0};
  constexpr std::array<double, 5> distances{-1.0, 0.0, 1.25, 2.0, 5.0};
  unsigned int bondRows = 0;
  unsigned int bondCalls = 0;
  unsigned int disagreements = 0;

  for (double r0 : r0Values) {
    for (double kb : kbValues) {
      for (double distance : distances) {
        const double first =
            ForceFields::MMFF::Utils::calcBondStretchEnergy(r0, kb, distance);
        ++bondCalls;
        const double second =
            ForceFields::MMFF::Utils::calcBondStretchEnergy(r0, kb, distance);
        ++bondCalls;
        if (bits(first) != bits(second)) {
          std::fprintf(stderr, "bond repeat mismatch row=%u\n", bondRows + 1);
          ++disagreements;
        }
        emitBond(++bondRows, r0, kb, distance, first);
      }
    }
  }

  constexpr std::array<double, 3> theta0Values{0.0, 109.47, 180.0};
  constexpr std::array<double, 3> kaValues{-0.0, 0.5, 4.0};
  constexpr std::array<bool, 2> linearValues{false, true};
  constexpr std::array<double, 5> cosThetaValues{-1.0, -0.5, 0.0, 0.5, 1.0};
  unsigned int angleRows = 0;
  unsigned int angleCalls = 0;

  for (double theta0 : theta0Values) {
    for (double ka : kaValues) {
      for (bool isLinear : linearValues) {
        for (double cosTheta : cosThetaValues) {
          const double first = ForceFields::MMFF::Utils::calcAngleBendEnergy(
              theta0, ka, isLinear, cosTheta);
          ++angleCalls;
          const double second = ForceFields::MMFF::Utils::calcAngleBendEnergy(
              theta0, ka, isLinear, cosTheta);
          ++angleCalls;
          if (bits(first) != bits(second)) {
            std::fprintf(stderr, "angle repeat mismatch row=%u\n",
                         angleRows + 1);
            ++disagreements;
          }
          emitAngle(++angleRows, theta0, ka, isLinear, cosTheta, first);
        }
      }
    }
  }

  constexpr std::array<std::uint64_t, 2> outOfDomainCosineBits{
      0xbff0000000000001ULL, 0x3ff0000000000001ULL};
  const double supplementTheta0 = 109.47;
  const double supplementKa = 0.5;
  for (std::uint64_t cosineBits : outOfDomainCosineBits) {
    const double cosTheta = std::bit_cast<double>(cosineBits);
    for (bool isLinear : linearValues) {
      const double first = ForceFields::MMFF::Utils::calcAngleBendEnergy(
          supplementTheta0, supplementKa, isLinear, cosTheta);
      ++angleCalls;
      emitAngle(++angleRows, supplementTheta0, supplementKa, isLinear,
                cosTheta, first);
    }
  }

  std::printf("COUNTS\tB\t%u\t%u\tA\t%u\t%u\n", bondRows, bondCalls,
              angleRows, angleCalls);
  if (bondRows != 45 || bondCalls != 90 || angleRows != 94 ||
      angleCalls != 184 || disagreements != 0) {
    std::fprintf(stderr,
                 "count/repeat failure B=%u/%u A=%u/%u repeat=%u\n",
                 bondRows, bondCalls, angleRows, angleCalls, disagreements);
    return 1;
  }
  return 0;
}
