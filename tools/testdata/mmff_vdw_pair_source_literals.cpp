#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

using std::exp;
using std::sqrt;

struct MMFFVdWCollection {
  double B;
  double Beta;
  double DARAD;
  double DAEPS;
};

// Keep only fields read by the three copied source functions.
struct MMFFVdW {
  double alpha_i;
  double N_i;
  double G_i;
  double R_star;
  std::uint8_t DA;
};

double calcUnscaledVdWMinimum(const MMFFVdWCollection *mmffVdW,
                              const MMFFVdW *mmffVdWParamsIAtom,
                              const MMFFVdW *mmffVdWParamsJAtom) {
  double gamma_ij = (mmffVdWParamsIAtom->R_star - mmffVdWParamsJAtom->R_star) /
                    (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star);

  return (0.5 * (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star) *
          (1.0 +
           (((mmffVdWParamsIAtom->DA == 'D') || (mmffVdWParamsJAtom->DA == 'D'))
                ? 0.0
                : mmffVdW->B *
                      (1.0 - exp(-(mmffVdW->Beta) * gamma_ij * gamma_ij)))));
}

double calcUnscaledVdWWellDepth(double R_star_ij,
                                const MMFFVdW *mmffVdWParamsIAtom,
                                const MMFFVdW *mmffVdWParamsJAtom) {
  double R_star_ij2 = R_star_ij * R_star_ij;
  double const c4 = 181.16;

  return (c4 * mmffVdWParamsIAtom->G_i * mmffVdWParamsJAtom->G_i *
          mmffVdWParamsIAtom->alpha_i * mmffVdWParamsJAtom->alpha_i /
          ((sqrt(mmffVdWParamsIAtom->alpha_i / mmffVdWParamsIAtom->N_i) +
            sqrt(mmffVdWParamsJAtom->alpha_i / mmffVdWParamsJAtom->N_i)) *
           R_star_ij2 * R_star_ij2 * R_star_ij2));
}

void scaleVdWParams(double &R_star_ij, double &wellDepth,
                    const MMFFVdWCollection *mmffVdW,
                    const MMFFVdW *mmffVdWParamsIAtom,
                    const MMFFVdW *mmffVdWParamsJAtom) {
  if (((mmffVdWParamsIAtom->DA == 'D') && (mmffVdWParamsJAtom->DA == 'A')) ||
      ((mmffVdWParamsIAtom->DA == 'A') && (mmffVdWParamsJAtom->DA == 'D'))) {
    R_star_ij *= mmffVdW->DARAD;
    wellDepth *= mmffVdW->DAEPS;
  }
}

struct Profile {
  MMFFVdWCollection collection;
  double i_alpha;
  double i_N;
  double i_A;
  double i_G;
  double i_R;
  double j_alpha;
  double j_N;
  double j_A;
  double j_G;
  double j_R;
  double depth_radius;
  double scale_radius;
  double scale_depth;
};

static unsigned long long bits(double value) {
  return static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(value));
}

static MMFFVdW row(double alpha, double n, double g, double r, std::uint8_t da) {
  return MMFFVdW{alpha, n, g, r, da};
}

static void emitProfile(unsigned int index, const Profile &p) {
  std::printf(
      "PROFILE\t%u\tB=%016llx\tBeta=%016llx\tDARAD=%016llx\tDAEPS=%016llx"
      "\ti(alpha,N,A,G,R)=%016llx,%016llx,%016llx,%016llx,%016llx"
      "\tj(alpha,N,A,G,R)=%016llx,%016llx,%016llx,%016llx,%016llx"
      "\tdepthR=%016llx\tscaleR=%016llx\tscaleDepth=%016llx\n",
      index, bits(p.collection.B), bits(p.collection.Beta),
      bits(p.collection.DARAD), bits(p.collection.DAEPS), bits(p.i_alpha),
      bits(p.i_N), bits(p.i_A), bits(p.i_G), bits(p.i_R), bits(p.j_alpha),
      bits(p.j_N), bits(p.j_A), bits(p.j_G), bits(p.j_R),
      bits(p.depth_radius), bits(p.scale_radius), bits(p.scale_depth));
}

static void emitBase(unsigned int index, unsigned int profileIndex,
                     const Profile &p, std::uint8_t iDa, std::uint8_t jDa) {
  const MMFFVdW i = row(p.i_alpha, p.i_N, p.i_G, p.i_R, iDa);
  const MMFFVdW j = row(p.j_alpha, p.j_N, p.j_G, p.j_R, jDa);
  const double minimum = calcUnscaledVdWMinimum(&p.collection, &i, &j);
  const double depth =
      calcUnscaledVdWWellDepth(p.depth_radius, &i, &j);
  double scaledRadius = p.scale_radius;
  double scaledDepth = p.scale_depth;
  scaleVdWParams(scaledRadius, scaledDepth, &p.collection, &i, &j);
  std::printf(
      "BASE\t%02u\tprofile=%u\tDA=%c,%c\tmin=%016llx"
      "\tdepthIn=%016llx\tdepth=%016llx"
      "\tscaleIn=%016llx,%016llx\tscaleOut=%016llx,%016llx\n",
      index, profileIndex, iDa, jDa, bits(minimum), bits(p.depth_radius),
      bits(depth), bits(p.scale_radius), bits(p.scale_depth), bits(scaledRadius),
      bits(scaledDepth));
}

static void emitScaleControl(unsigned int index, const Profile &p,
                             std::uint8_t iDa, std::uint8_t jDa,
                             double radius, double depth) {
  const MMFFVdW i = row(p.i_alpha, p.i_N, p.i_G, p.i_R, iDa);
  const MMFFVdW j = row(p.j_alpha, p.j_N, p.j_G, p.j_R, jDa);
  const double inputRadius = radius;
  const double inputDepth = depth;
  scaleVdWParams(radius, depth, &p.collection, &i, &j);
  std::printf(
      "SCALE\t%02u\tDA=%c,%c\tin=%016llx,%016llx"
      "\tout=%016llx,%016llx\n",
      index, iDa, jDa, bits(inputRadius), bits(inputDepth), bits(radius),
      bits(depth));
}

int main() {
  constexpr std::array<Profile, 2> profiles{{
      {{0.2, 12.0, 0.8, 0.5}, 4.0, 2.0, 3.0, 1.25, 3.0, 9.0, 3.0,
       5.0, 0.75, 5.0, 2.75, 1.75, 0.625},
      {{-0.0, -0.0, -0.0, -2.0}, 16.0, 4.0, 4.0, 2.0, 4.0, 1.0, 1.0,
       4.0, 0.5, 4.0, 3.5, 2.25, -1.5},
  }};
  constexpr std::array<std::uint8_t, 4> daValues{'D', 'A', '-', 'X'};
  unsigned int index = 0;
  for (unsigned int profileIndex = 0; profileIndex < profiles.size();
       ++profileIndex) {
    const Profile &p = profiles[profileIndex];
    emitProfile(profileIndex, p);
    for (const std::uint8_t iDa : daValues) {
      for (const std::uint8_t jDa : daValues) {
        emitBase(index++, profileIndex, p, iDa, jDa);
      }
    }
  }

  constexpr std::array<std::array<std::uint8_t, 2>, 4> scaleDaPairs{{
      {{'D', 'A'}}, {{'A', 'D'}}, {{'D', 'D'}}, {{'-', 'A'}},
  }};
  const Profile &p = profiles[0];
  index = 0;
  for (const auto &pair : scaleDaPairs) {
    emitScaleControl(index++, p, pair[0], pair[1], -0.0, +0.0);
    emitScaleControl(index++, p, pair[0], pair[1],
                     std::numeric_limits<double>::infinity(),
                     -std::numeric_limits<double>::infinity());
  }
}
