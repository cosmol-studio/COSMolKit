#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>

using std::pow;

struct Cell {
  const char *family;
  unsigned int index;
  double alpha;
  double power;
  double a;
};

static void emit(const Cell &cell) {
  // Pinned RDKit Params.cpp: MMFFVdWCollection constructor computes
  // mmffVdWObj.A_i * pow(mmffVdWObj.alpha_i, this->power).
  const double r_star = cell.a * pow(cell.alpha, cell.power);
  std::printf(
      "%s\t%u\t%a\t%016llx\t%a\t%016llx\t%a\t%016llx\t%016llx\n",
      cell.family, cell.index, cell.alpha,
      static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(cell.alpha)),
      cell.power,
      static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(cell.power)),
      cell.a,
      static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(cell.a)),
      static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(r_star)));
}

int main() {
  constexpr std::array<double, 3> alphas{1.0, 16.0, 81.0};
  constexpr std::array<double, 4> powers{0.0, 0.25, 0.5, 1.0};
  unsigned int index = 0;
  for (const double alpha : alphas) {
    for (const double power : powers) {
      emit(Cell{"cartesian", index++, alpha, power, 2.0});
    }
  }

  constexpr std::array<Cell, 4> source_sentinels{{
      {"default", 1, 1.050, 0.25, 3.890},
      {"default", 21, 0.150, 0.25, 4.200},
      {"default", 82, 0.950, 0.25, 3.890},
      {"default", 99, 0.35, 0.25, 4.0},
  }};
  for (const Cell &cell : source_sentinels) {
    emit(cell);
  }
}
