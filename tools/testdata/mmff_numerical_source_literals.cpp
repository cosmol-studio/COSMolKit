// Source-only pinned RDKit numerical reference for MMFF-NUMERICAL1-62.
// Pin: 351f8f378f8ad6bbd517980c38896e66bf907af8c.
// The selected bodies below are copied from the named source functions. This
// program has no COSMolKit link or dependency and emits fixed source products.
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <tuple>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RDGeom {

// RDKit Code/Geometry/point.h: Point3D coordinate fields and constructor.
struct Point3D {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  constexpr Point3D() = default;
  constexpr Point3D(double xv, double yv, double zv) : x(xv), y(yv), z(zv) {}

  // RDKit point.h:132-137, operator/=; exact component order and expressions.
  constexpr Point3D &operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
  }

  // RDKit point.h:139-145, unary operator-.
  constexpr Point3D operator-() const {
    Point3D res(x, y, z);
    res.x *= -1.0;
    res.y *= -1.0;
    res.z *= -1.0;
    return res;
  }

  // RDKit point.h:158-161, length.
  double length() const {
    double res = x * x + y * y + z * z;
    return sqrt(res);
  }

  // RDKit point.h:169-172, dotProduct.
  constexpr double dotProduct(const Point3D &other) const {
    double res = x * (other.x) + y * (other.y) + z * (other.z);
    return res;
  }

  // RDKit point.h:228-234, crossProduct.
  constexpr Point3D crossProduct(const Point3D &other) const {
    Point3D res;
    res.x = y * (other.z) - z * (other.y);
    res.y = -x * (other.z) + z * (other.x);
    res.z = x * (other.y) - y * (other.x);
    return res;
  }
};

// RDKit Code/Geometry/point.cpp:64-70, binary operator-.
Point3D operator-(const Point3D &p1, const Point3D &p2) {
  Point3D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  res.z = p1.z - p2.z;
  return res;
}

}  // namespace RDGeom

// The only source PRECONDITION admissions exercised are non-null. Invalid
// native pointers are excluded; this driver guard makes accidental misuse fail.
#define PRECONDITION(condition, message) \
  do {                                  \
    if (!(condition)) {                 \
      std::fprintf(stderr, "invalid excluded source precondition: %s\n", message); \
      std::abort();                     \
    }                                   \
  } while (false)

namespace ForceFields {
namespace MMFF {

constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180.0 / M_PI;
constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;

inline bool isDoubleZero(const double x) {
  return ((x < 1.0e-10) && (x > -1.0e-10));
}

inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }

struct MMFFBond {
  double kb;
  double r0;
};
struct MMFFAngle {
  double ka;
  double theta0;
};
struct MMFFStbn {
  double kbaIJK;
  double kbaKJI;
};
struct MMFFTor {
  double V1;
  double V2;
  double V3;
};
struct MMFFOop {
  double koop;
};

namespace Utils {

// BEGIN RDKit source BondStretch.cpp:20-24 calcBondRestLength.
double calcBondRestLength(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->r0;
}
// END RDKit source calcBondRestLength.

// BEGIN RDKit source BondStretch.cpp:26-30 calcBondForceConstant.
double calcBondForceConstant(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->kb;
}
// END RDKit source calcBondForceConstant.

// BEGIN RDKit source AngleBend.cpp:21-25 calcAngleRestValue.
double calcAngleRestValue(const MMFFAngle *mmffAngleParams) {
  PRECONDITION(mmffAngleParams, "angle parameters not found");

  return mmffAngleParams->theta0;
}
// END RDKit source calcAngleRestValue.

// BEGIN RDKit source AngleBend.cpp:27-35 calcCosTheta.
double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2, RDGeom::Point3D p3,
                    double dist1, double dist2) {
  RDGeom::Point3D p12 = p1 - p2;
  RDGeom::Point3D p32 = p3 - p2;
  double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
  clipToOne(cosTheta);

  return cosTheta;
}
// END RDKit source calcCosTheta.

// BEGIN RDKit source AngleBend.cpp:37-41 calcAngleForceConstant.
double calcAngleForceConstant(const MMFFAngle *mmffAngleParams) {
  PRECONDITION(mmffAngleParams, "angle parameters not found");

  return mmffAngleParams->ka;
}
// END RDKit source calcAngleForceConstant.

// BEGIN RDKit source StretchBend.cpp:22-27 calcStbnForceConstants.
std::pair<double, double> calcStbnForceConstants(
    const MMFFStbn *mmffStbnParams) {
  PRECONDITION(mmffStbnParams, "stretch-bend parameters not found");

  return std::make_pair(mmffStbnParams->kbaIJK, mmffStbnParams->kbaKJI);
}
// END RDKit source calcStbnForceConstants.

// BEGIN RDKit source TorsionAngle.cpp:40-44 calcTorsionForceConstant.
std::tuple<double, double, double> calcTorsionForceConstant(
    const MMFFTor *mmffTorParams) {
  return std::make_tuple(mmffTorParams->V1, mmffTorParams->V2,
                         mmffTorParams->V3);
}
// END RDKit source calcTorsionForceConstant.

// BEGIN RDKit source OopBend.cpp:36-40 calcOopBendForceConstant.
double calcOopBendForceConstant(const MMFFOop *mmffOopParams) {
  PRECONDITION(mmffOopParams, "no OOP parameters");

  return mmffOopParams->koop;
}
// END RDKit source calcOopBendForceConstant.

// BEGIN RDKit source StretchBend.cpp:29-37 calcStretchBendEnergy.
std::pair<double, double> calcStretchBendEnergy(
    const double deltaDist1, const double deltaDist2, const double deltaTheta,
    const std::pair<double, double> forceConstants) {
  double factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * deltaTheta;

  return std::make_pair(factor * forceConstants.first * deltaDist1,
                        factor * forceConstants.second * deltaDist2);
}
// END RDKit source calcStretchBendEnergy.

// BEGIN RDKit source TorsionAngle.cpp:46-53 calcTorsionEnergy.
double calcTorsionEnergy(const double V1, const double V2, const double V3,
                         const double cosPhi) {
  double cos2Phi = 2.0 * cosPhi * cosPhi - 1.0;
  double cos3Phi = cosPhi * (2.0 * cos2Phi - 1.0);

  return (0.5 *
          (V1 * (1.0 + cosPhi) + V2 * (1.0 - cos2Phi) + V3 * (1.0 + cos3Phi)));
}
// END RDKit source calcTorsionEnergy.

// BEGIN RDKit source OopBend.cpp:42-45 calcOopBendEnergy.
double calcOopBendEnergy(const double chi, const double koop) {
  double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
  return (0.5 * c2 * koop * chi * chi);
}
// END RDKit source calcOopBendEnergy.

// BEGIN RDKit source TorsionAngle.cpp:19-38 calcTorsionCosPhi.
double calcTorsionCosPhi(const RDGeom::Point3D &iPoint,
                         const RDGeom::Point3D &jPoint,
                         const RDGeom::Point3D &kPoint,
                         const RDGeom::Point3D &lPoint) {
  RDGeom::Point3D r1 = iPoint - jPoint;
  RDGeom::Point3D r2 = kPoint - jPoint;
  RDGeom::Point3D r3 = jPoint - kPoint;
  RDGeom::Point3D r4 = lPoint - kPoint;
  RDGeom::Point3D t1 = r1.crossProduct(r2);
  RDGeom::Point3D t2 = r3.crossProduct(r4);
  auto t1_len = t1.length();
  auto t2_len = t2.length();
  if (isDoubleZero(t1_len) || isDoubleZero(t2_len)) {
    return 0.0;
  }
  double cosPhi = t1.dotProduct(t2) / (t1_len * t2_len);
  clipToOne(cosPhi);

  return cosPhi;
}
// END RDKit source calcTorsionCosPhi.

// BEGIN RDKit source OopBend.cpp:18-34 calcOopChi.
double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
                  const RDGeom::Point3D &kPoint,
                  const RDGeom::Point3D &lPoint) {
  RDGeom::Point3D rJI = iPoint - jPoint;
  RDGeom::Point3D rJK = kPoint - jPoint;
  RDGeom::Point3D rJL = lPoint - jPoint;
  rJI /= rJI.length();
  rJK /= rJK.length();
  rJL /= rJL.length();

  RDGeom::Point3D n = rJI.crossProduct(rJK);
  n /= n.length();
  double sinChi = n.dotProduct(rJL);
  clipToOne(sinChi);

  return RAD2DEG * asin(sinChi);
}
// END RDKit source calcOopChi.

// BEGIN RDKit source AngleBend.cpp:59-86 calcAngleBendGrad.
void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
                       double &dE_dTheta, double &cosTheta, double &sinTheta) {
  // -------
  // dTheta/dx is trickier:
  double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
                       1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
                       1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
                       1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
                       1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
                       1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};

  g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
  g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
  g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);

  g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
  g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
  g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);

  g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
  g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
  g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
}
// END RDKit source calcAngleBendGrad.

// BEGIN RDKit source TorsionAngle.cpp:55-93 calcTorsionGrad.
void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
                     double **g, double &sinTerm, double &cosPhi) {
  // -------
  // dTheta/dx is trickier:
  double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
                       1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
                       1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
                       1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
                       1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
                       1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};

  g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
  g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
  g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);

  g[1][0] += sinTerm *
             (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
              dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
  g[1][1] += sinTerm *
             (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
              dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
  g[1][2] += sinTerm *
             (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
              dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));

  g[2][0] += sinTerm *
             (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
              dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
  g[2][1] += sinTerm *
             (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
              dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
  g[2][2] += sinTerm *
             (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
              dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));

  g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
  g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
  g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
}
// END RDKit source calcTorsionGrad.

}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields

namespace {

using ForceFields::MMFF::MMFFAngle;
using ForceFields::MMFF::MMFFBond;
using ForceFields::MMFF::MMFFOop;
using ForceFields::MMFF::MMFFStbn;
using ForceFields::MMFF::MMFFTor;
using RDGeom::Point3D;

constexpr std::array<std::size_t, 9> EXPECTED_CELLS = {75, 108, 48, 18, 24,
                                                       20, 16, 96, 192};
constexpr std::array<std::size_t, 9> EXPECTED_CALLS = {150, 216, 96, 36, 48,
                                                       40, 32, 192, 384};
std::array<std::size_t, 9> cells{};
std::array<std::size_t, 9> calls{};
std::size_t repeat_mismatches = 0;

std::uint64_t bits(double value) {
  return std::bit_cast<std::uint64_t>(value);
}

double from_bits(std::uint64_t value) {
  return std::bit_cast<double>(value);
}

bool same_bits(double left, double right) { return bits(left) == bits(right); }

bool same_bits(const std::pair<double, double> &left,
               const std::pair<double, double> &right) {
  return same_bits(left.first, right.first) &&
         same_bits(left.second, right.second);
}

bool same_bits(const std::tuple<double, double, double> &left,
               const std::tuple<double, double, double> &right) {
  return same_bits(std::get<0>(left), std::get<0>(right)) &&
         same_bits(std::get<1>(left), std::get<1>(right)) &&
         same_bits(std::get<2>(left), std::get<2>(right));
}

void start_row(unsigned int unit, const char *label) {
  std::printf("U%u\t%s", unit, label);
}

void print_value(double value) {
  std::printf("\t%016llx", static_cast<unsigned long long>(bits(value)));
}

void print_result(unsigned int unit, const char *label, double result) {
  start_row(unit, label);
  print_value(result);
  std::printf("\n");
}

void print_result(unsigned int unit, const char *label,
                  const std::pair<double, double> &result) {
  start_row(unit, label);
  print_value(result.first);
  print_value(result.second);
  std::printf("\n");
}

void print_result(unsigned int unit, const char *label,
                  const std::tuple<double, double, double> &result) {
  start_row(unit, label);
  print_value(std::get<0>(result));
  print_value(std::get<1>(result));
  print_value(std::get<2>(result));
  std::printf("\n");
}

template <class Function>
auto twice(unsigned int unit, const char *label, Function &&function) {
  auto first = function();
  ++calls[unit - 1];
  auto second = function();
  ++calls[unit - 1];
  ++cells[unit - 1];
  if (!same_bits(first, second)) {
    ++repeat_mismatches;
    std::fprintf(stderr, "repeat mismatch U%u %s\n", unit, label);
  }
  print_result(unit, label, first);
  return first;
}

using PointBits = std::array<std::array<std::uint64_t, 3>, 4>;

constexpr std::uint64_t Q_BITS[10][4][3] = {
    {{0x3ff0000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3ff0000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3ff0000000000000ULL, 0x3ff0000000000000ULL}},
    {{0x3ff0000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0xbff0000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0xbff0000000000000ULL, 0x3ff0000000000000ULL, 0x0000000000000000ULL}},
    {{0x3ff0000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x3ff0000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x0000000000000000ULL}},
    {{0x3fd0000000000000ULL, 0xbfe0000000000000ULL, 0x3ff0000000000000ULL},
     {0xbff4000000000000ULL, 0x3fe8000000000000ULL, 0x4000000000000000ULL},
     {0x3fe0000000000000ULL, 0x3ff8000000000000ULL, 0xbfe8000000000000ULL},
     {0x4000000000000000ULL, 0xbfd0000000000000ULL, 0x3fc0000000000000ULL}},
    {{0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL}},
    {{0x4010000000000000ULL, 0xc000000000000000ULL, 0x3fe0000000000000ULL},
     {0x4008000000000000ULL, 0xc000000000000000ULL, 0x3fe0000000000000ULL},
     {0x4008000000000000ULL, 0xbff0000000000000ULL, 0x3fe0000000000000ULL},
     {0x4008000000000000ULL, 0xbff0000000000000ULL, 0x3ff8000000000000ULL}},
    {{0x3fa0000000000000ULL, 0xbfb0000000000000ULL, 0x3fc0000000000000ULL},
     {0xbfc4000000000000ULL, 0x3fb8000000000000ULL, 0x3fd0000000000000ULL},
     {0x3fb0000000000000ULL, 0x3fc8000000000000ULL, 0xbfb8000000000000ULL},
     {0x3fd0000000000000ULL, 0xbfa0000000000000ULL, 0x3f90000000000000ULL}},
    {{0x3ff0000000000000ULL, 0x8000000000000000ULL, 0x8000000000000000ULL},
     {0x8000000000000000ULL, 0x8000000000000000ULL, 0x8000000000000000ULL},
     {0x8000000000000000ULL, 0x3ff0000000000000ULL, 0x8000000000000000ULL},
     {0x8000000000000000ULL, 0x3ff0000000000000ULL, 0x3ff0000000000000ULL}},
    {{0x3eb0c6f7a0b5ed8dULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3eb0c6f7a0b5ed8dULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3eb0c6f7a0b5ed8dULL, 0x3eb0c6f7a0b5ed8dULL}},
    {{0x3f1a36e2eb1c432dULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3f1a36e2eb1c432dULL, 0x0000000000000000ULL},
     {0x0000000000000000ULL, 0x3f1a36e2eb1c432dULL, 0x3f1a36e2eb1c432dULL}},
};

std::array<Point3D, 4> points(std::size_t set) {
  std::array<Point3D, 4> result{};
  for (std::size_t row = 0; row < result.size(); ++row) {
    result[row] = Point3D(from_bits(Q_BITS[set][row][0]),
                          from_bits(Q_BITS[set][row][1]),
                          from_bits(Q_BITS[set][row][2]));
  }
  return result;
}

std::array<Point3D, 4> reversed(std::array<Point3D, 4> value) {
  std::reverse(value.begin(), value.end());
  return value;
}

using Gradient = std::array<std::array<double, 3>, 7>;
constexpr std::uint64_t S_BITS[7][3] = {
    {0x3fd0000000000000ULL, 0xbfe0000000000000ULL, 0x8000000000000000ULL},
    {0x3fe0000000000000ULL, 0xbff0000000000000ULL, 0x8000000000000000ULL},
    {0x3fe8000000000000ULL, 0xbff8000000000000ULL, 0x8000000000000000ULL},
    {0x3ff0000000000000ULL, 0xc000000000000000ULL, 0x8000000000000000ULL},
    {0x3ff4000000000000ULL, 0xc004000000000000ULL, 0x8000000000000000ULL},
    {0x3ff8000000000000ULL, 0xc008000000000000ULL, 0x8000000000000000ULL},
    {0x3ffc000000000000ULL, 0xc00c000000000000ULL, 0x8000000000000000ULL},
};

Gradient baseline(bool nonzero) {
  Gradient result{};
  if (nonzero) {
    for (std::size_t row = 0; row < result.size(); ++row) {
      for (std::size_t column = 0; column < 3; ++column) {
        result[row][column] = from_bits(S_BITS[row][column]);
      }
    }
  }
  return result;
}

std::array<std::uint64_t, 21> gradient_bits(const Gradient &gradient) {
  std::array<std::uint64_t, 21> result{};
  std::size_t index = 0;
  for (const auto &row : gradient) {
    for (double value : row) {
      result[index++] = bits(value);
    }
  }
  return result;
}

void print_gradient(unsigned int unit, const char *label,
                    const std::array<std::uint64_t, 21> &result) {
  start_row(unit, label);
  for (std::uint64_t value : result) {
    std::printf("\t%016llx", static_cast<unsigned long long>(value));
  }
  std::printf("\n");
}

void run_angle_gradient(const char *label, const std::array<Point3D, 2> &r_input,
                        const std::array<double, 2> &dist_input,
                        double d_e_input, double cos_input, double sin_input,
                        const std::array<std::size_t, 3> &indices,
                        bool nonzero) {
  std::array<std::uint64_t, 21> first{};
  for (unsigned int repeat = 0; repeat < 2; ++repeat) {
    Gradient gradient = baseline(nonzero);
    Point3D r[2] = {r_input[0], r_input[1]};
    double dist[2] = {dist_input[0], dist_input[1]};
    double d_e_d_theta = d_e_input;
    double cos_theta = cos_input;
    double sin_theta = sin_input;
    double *g[3] = {gradient[indices[0]].data(), gradient[indices[1]].data(),
                    gradient[indices[2]].data()};
    ForceFields::MMFF::Utils::calcAngleBendGrad(
        r, dist, g, d_e_d_theta, cos_theta, sin_theta);
    ++calls[7];
    const auto current = gradient_bits(gradient);
    if (repeat == 0) {
      first = current;
    } else if (current != first) {
      ++repeat_mismatches;
      std::fprintf(stderr, "repeat mismatch U8 %s\n", label);
    }
  }
  ++cells[7];
  print_gradient(8, label, first);
}

void run_torsion_gradient(const char *label,
                          const std::array<Point3D, 4> &r_input,
                          const std::array<Point3D, 2> &t_input,
                          const std::array<double, 2> &d_input,
                          double sin_input, double cos_input,
                          const std::array<std::size_t, 4> &indices,
                          bool nonzero) {
  std::array<std::uint64_t, 21> first{};
  for (unsigned int repeat = 0; repeat < 2; ++repeat) {
    Gradient gradient = baseline(nonzero);
    Point3D r[4] = {r_input[0], r_input[1], r_input[2], r_input[3]};
    Point3D t[2] = {t_input[0], t_input[1]};
    double d[2] = {d_input[0], d_input[1]};
    double sin_term = sin_input;
    double cos_phi = cos_input;
    double *g[4] = {gradient[indices[0]].data(), gradient[indices[1]].data(),
                    gradient[indices[2]].data(), gradient[indices[3]].data()};
    ForceFields::MMFF::Utils::calcTorsionGrad(r, t, d, g, sin_term, cos_phi);
    ++calls[8];
    const auto current = gradient_bits(gradient);
    if (repeat == 0) {
      first = current;
    } else if (current != first) {
      ++repeat_mismatches;
      std::fprintf(stderr, "repeat mismatch U9 %s\n", label);
    }
  }
  ++cells[8];
  print_gradient(9, label, first);
}

}  // namespace

int main() {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  static_assert(std::numeric_limits<double>::is_iec559);

  constexpr std::uint64_t V_BITS[] = {0x8000000000000000ULL,
                                      0x3fe0000000000000ULL,
                                      0x4010000000000000ULL};
  std::array<double, 3> v = {from_bits(V_BITS[0]), from_bits(V_BITS[1]),
                             from_bits(V_BITS[2])};
  char label[64];
  std::size_t ordinal = 0;

  for (std::size_t r0 = 0; r0 < 3; ++r0) {
    for (std::size_t kb = 0; kb < 3; ++kb) {
      MMFFBond params{v[kb], v[r0]};
      std::snprintf(label, sizeof(label), "bond_r0_%02zu", ordinal);
      twice(1, label, [&]() {
        return ForceFields::MMFF::Utils::calcBondRestLength(&params);
      });
      std::snprintf(label, sizeof(label), "bond_kb_%02zu", ordinal);
      twice(1, label, [&]() {
        return ForceFields::MMFF::Utils::calcBondForceConstant(&params);
      });
      ++ordinal;
    }
  }

  ordinal = 0;
  for (std::size_t theta0 = 0; theta0 < 3; ++theta0) {
    for (std::size_t ka = 0; ka < 3; ++ka) {
      MMFFAngle params{v[ka], v[theta0]};
      std::snprintf(label, sizeof(label), "angle_theta0_%02zu", ordinal);
      twice(1, label, [&]() {
        return ForceFields::MMFF::Utils::calcAngleRestValue(&params);
      });
      std::snprintf(label, sizeof(label), "angle_ka_%02zu", ordinal);
      twice(1, label, [&]() {
        return ForceFields::MMFF::Utils::calcAngleForceConstant(&params);
      });
      ++ordinal;
    }
  }

  ordinal = 0;
  for (double first : v) {
    for (double second : v) {
      MMFFStbn params{first, second};
      std::snprintf(label, sizeof(label), "stbn_%02zu", ordinal++);
      twice(1, label, [&]() {
        return ForceFields::MMFF::Utils::calcStbnForceConstants(&params);
      });
    }
  }

  ordinal = 0;
  for (double v1 : v) {
    for (double v2 : v) {
      for (double v3 : v) {
        MMFFTor params{v1, v2, v3};
        std::snprintf(label, sizeof(label), "tor_%02zu", ordinal++);
        twice(1, label, [&]() {
          return ForceFields::MMFF::Utils::calcTorsionForceConstant(&params);
        });
      }
    }
  }

  ordinal = 0;
  for (double koop : v) {
    MMFFOop params{koop};
    std::snprintf(label, sizeof(label), "oop_%02zu", ordinal++);
    twice(1, label, [&]() {
      return ForceFields::MMFF::Utils::calcOopBendForceConstant(&params);
    });
  }

  constexpr std::uint64_t U2_D1[] = {0xbfc0000000000000ULL,
                                     0x0000000000000000ULL,
                                     0x3fe8000000000000ULL};
  constexpr std::uint64_t U2_D2[] = {0xbfe0000000000000ULL,
                                     0x0000000000000000ULL,
                                     0x3ff4000000000000ULL};
  constexpr std::uint64_t U2_THETA[] = {0xc03e000000000000ULL,
                                        0x8000000000000000ULL,
                                        0x4046800000000000ULL};
  constexpr std::uint64_t U2_K1[] = {0x8000000000000000ULL,
                                     0x3fe0000000000000ULL};
  constexpr std::uint64_t U2_K2[] = {0xbfe8000000000000ULL,
                                     0x3ff4000000000000ULL};
  ordinal = 0;
  for (std::uint64_t d1 : U2_D1) {
    for (std::uint64_t d2 : U2_D2) {
      for (std::uint64_t theta : U2_THETA) {
        for (std::uint64_t k1 : U2_K1) {
          for (std::uint64_t k2 : U2_K2) {
            std::snprintf(label, sizeof(label), "stretch_bend_%03zu", ordinal++);
            twice(2, label, [&]() {
              return ForceFields::MMFF::Utils::calcStretchBendEnergy(
                  from_bits(d1), from_bits(d2), from_bits(theta),
                  {from_bits(k1), from_bits(k2)});
            });
          }
        }
      }
    }
  }

  constexpr std::uint64_t U3_V1[] = {0x8000000000000000ULL,
                                     0x3fe0000000000000ULL};
  constexpr std::uint64_t U3_V2[] = {0xbff4000000000000ULL,
                                     0x3fe8000000000000ULL};
  constexpr std::uint64_t U3_V3[] = {0x0000000000000000ULL,
                                     0x4000000000000000ULL};
  constexpr std::uint64_t U3_COS[] = {0xbff0000000000000ULL,
                                      0xbfe0000000000000ULL,
                                      0x8000000000000000ULL,
                                      0x3fe0000000000000ULL,
                                      0x3ff0000000000000ULL,
                                      0x3ff4000000000000ULL};
  ordinal = 0;
  for (std::uint64_t v1 : U3_V1) {
    for (std::uint64_t v2 : U3_V2) {
      for (std::uint64_t v3 : U3_V3) {
        for (std::uint64_t cos_phi : U3_COS) {
          std::snprintf(label, sizeof(label), "torsion_energy_%03zu", ordinal++);
          twice(3, label, [&]() {
            return ForceFields::MMFF::Utils::calcTorsionEnergy(
                from_bits(v1), from_bits(v2), from_bits(v3), from_bits(cos_phi));
          });
        }
      }
    }
  }

  constexpr std::uint64_t U4_CHI[] = {0xc05e000000000000ULL,
                                      0xc03e000000000000ULL,
                                      0x8000000000000000ULL,
                                      0x0000000000000000ULL,
                                      0x402e000000000000ULL,
                                      0x4056800000000000ULL};
  constexpr std::uint64_t U4_KOOP[] = {0x8000000000000000ULL,
                                       0x3fe0000000000000ULL,
                                       0x4000000000000000ULL};
  ordinal = 0;
  for (std::uint64_t chi : U4_CHI) {
    for (std::uint64_t koop : U4_KOOP) {
      std::snprintf(label, sizeof(label), "oop_energy_%02zu", ordinal++);
      twice(4, label, [&]() {
        return ForceFields::MMFF::Utils::calcOopBendEnergy(from_bits(chi),
                                                            from_bits(koop));
      });
    }
  }

  constexpr std::uint64_t U5_DIST[3][2] = {
      {0x3ff0000000000000ULL, 0x3ff0000000000000ULL},
      {0x3fe0000000000000ULL, 0x4000000000000000ULL},
      {0x0000000000000000ULL, 0x3ff0000000000000ULL}};
  for (std::size_t q = 0; q < 8; ++q) {
    const auto pts = points(q);
    for (std::size_t d = 0; d < 3; ++d) {
      std::snprintf(label, sizeof(label), "cos_theta_q%zu_d%zu", q, d);
      twice(5, label, [&]() {
        return ForceFields::MMFF::Utils::calcCosTheta(
            pts[0], pts[1], pts[2], from_bits(U5_DIST[d][0]),
            from_bits(U5_DIST[d][1]));
      });
    }
  }

  for (std::size_t q = 0; q < 10; ++q) {
    const auto pts = points(q);
    for (unsigned int order = 0; order < 2; ++order) {
      const auto ordered = order == 0 ? pts : reversed(pts);
      std::snprintf(label, sizeof(label), "cos_phi_q%zu_order%u", q, order);
      twice(6, label, [&]() {
        return ForceFields::MMFF::Utils::calcTorsionCosPhi(
            ordered[0], ordered[1], ordered[2], ordered[3]);
      });
    }
  }

  for (std::size_t q = 0; q < 8; ++q) {
    const auto pts = points(q);
    for (unsigned int order = 0; order < 2; ++order) {
      const auto ordered = order == 0 ? pts : reversed(pts);
      std::snprintf(label, sizeof(label), "oop_chi_q%zu_order%u", q, order);
      twice(7, label, [&]() {
        return ForceFields::MMFF::Utils::calcOopChi(
            ordered[0], ordered[1], ordered[2], ordered[3]);
      });
    }
  }

  constexpr std::uint64_t U8_DIST[2][2] = {
      {0x3ff0000000000000ULL, 0x4000000000000000ULL},
      {0x3fe0000000000000ULL, 0x3ff4000000000000ULL}};
  constexpr std::uint64_t U8_DE[2] = {0x8000000000000000ULL,
                                       0x3ff4000000000000ULL};
  constexpr std::uint64_t U8_CS[2][2] = {
      {0xbfe0000000000000ULL, 0xbfe999999999999aULL},
      {0x3fe8000000000000ULL, 0x3fe3333333333333ULL}};
  constexpr std::size_t U8_INDICES[3][3] = {{1, 3, 5}, {2, 2, 2}, {1, 3, 1}};
  ordinal = 0;
  for (std::size_t r_set = 0; r_set < 2; ++r_set) {
    const auto pts = points(r_set == 0 ? 0 : 3);
    const std::array<Point3D, 2> r_input = {pts[0], pts[1]};
    for (std::size_t d = 0; d < 2; ++d) {
      for (std::uint64_t d_e : U8_DE) {
        for (std::size_t cs = 0; cs < 2; ++cs) {
          for (std::size_t ix = 0; ix < 3; ++ix) {
            for (unsigned int nonzero = 0; nonzero < 2; ++nonzero) {
              std::snprintf(label, sizeof(label), "angle_grad_%03zu", ordinal++);
              run_angle_gradient(
                  label, r_input,
                  {from_bits(U8_DIST[d][0]), from_bits(U8_DIST[d][1])},
                  from_bits(d_e), from_bits(U8_CS[cs][0]),
                  from_bits(U8_CS[cs][1]),
                  {U8_INDICES[ix][0], U8_INDICES[ix][1], U8_INDICES[ix][2]},
                  nonzero != 0);
            }
          }
        }
      }
    }
  }

  const std::array<Point3D, 2> U9_T[2] = {
      {Point3D(from_bits(0x3ff0000000000000ULL), 0.0, 0.0),
       Point3D(0.0, 1.0, 0.0)},
      {points(3)[0], points(3)[1]}};
  constexpr std::uint64_t U9_D[2][2] = {
      {0x3ff0000000000000ULL, 0x4000000000000000ULL},
      {0x3fe0000000000000ULL, 0x3ff4000000000000ULL}};
  constexpr std::uint64_t U9_SIN[2] = {0x8000000000000000ULL,
                                       0x3ff4000000000000ULL};
  constexpr std::uint64_t U9_COS[2] = {0xbfe0000000000000ULL,
                                       0x3fe8000000000000ULL};
  constexpr std::size_t U9_INDICES[3][4] = {
      {1, 3, 5, 6}, {2, 2, 2, 2}, {1, 3, 1, 3}};
  ordinal = 0;
  for (std::size_t r_set = 0; r_set < 2; ++r_set) {
    const auto r_input = points(r_set == 0 ? 0 : 3);
    for (std::size_t t_set = 0; t_set < 2; ++t_set) {
      for (std::size_t d = 0; d < 2; ++d) {
        for (std::uint64_t sin_term : U9_SIN) {
          for (std::uint64_t cos_phi : U9_COS) {
            for (std::size_t ix = 0; ix < 3; ++ix) {
              for (unsigned int nonzero = 0; nonzero < 2; ++nonzero) {
                std::snprintf(label, sizeof(label), "torsion_grad_%03zu", ordinal++);
                run_torsion_gradient(
                    label, r_input, U9_T[t_set],
                    {from_bits(U9_D[d][0]), from_bits(U9_D[d][1])},
                    from_bits(sin_term), from_bits(cos_phi),
                    {U9_INDICES[ix][0], U9_INDICES[ix][1], U9_INDICES[ix][2],
                     U9_INDICES[ix][3]},
                    nonzero != 0);
              }
            }
          }
        }
      }
    }
  }

  std::size_t total_cells = 0;
  std::size_t total_calls = 0;
  bool counts_match = true;
  for (std::size_t unit = 0; unit < cells.size(); ++unit) {
    total_cells += cells[unit];
    total_calls += calls[unit];
    if (cells[unit] != EXPECTED_CELLS[unit] ||
        calls[unit] != EXPECTED_CALLS[unit]) {
      counts_match = false;
      std::fprintf(stderr, "count mismatch U%zu cells=%zu calls=%zu\n", unit + 1,
                   cells[unit], calls[unit]);
    }
  }
  std::printf("#TOTAL\t%zu\t%zu\n", total_cells, total_calls);
  std::printf("#CONSTANTS\t%016llx\t%016llx\t%016llx\t%016llx\n",
              static_cast<unsigned long long>(bits(M_PI)),
              static_cast<unsigned long long>(bits(ForceFields::MMFF::DEG2RAD)),
              static_cast<unsigned long long>(bits(ForceFields::MMFF::RAD2DEG)),
              static_cast<unsigned long long>(
                  bits(ForceFields::MMFF::MDYNE_A_TO_KCAL_MOL)));
  if (!counts_match || total_cells != 597 || total_calls != 1194 ||
      repeat_mismatches != 0) {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
