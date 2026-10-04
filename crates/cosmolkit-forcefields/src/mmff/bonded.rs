pub(super) fn calc_bond_stretch_energy(r0: f64, kb: f64, distance: f64) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/BondStretch.cpp:32-42:
    // RDKit❗✔️: double calcBondStretchEnergy(const double r0, const double kb,
    // RDKit❗✔️:                              const double distance) {
    // RDKit❗✔️:   double distTerm = distance - r0;
    // RDKit❗✔️:   double distTerm2 = distTerm * distTerm;
    // RDKit❗✔️:   double const c1 = MDYNE_A_TO_KCAL_MOL;
    // RDKit❗✔️:   double const cs = -2.0;
    // RDKit❗✔️:   double const c3 = 7.0 / 12.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   return (0.5 * c1 * kb * distTerm2 *
    // RDKit❗✔️:           (1.0 + cs * distTerm + c3 * cs * cs * distTerm2));
    // RDKit❗✔️: }
    // RDKit❗✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
    // Behavior review: preserve the source subtraction, square, local constants,
    // left-associated product, and polynomial evaluation for every f64 input.
    // Complexity review: the owner performs a fixed number of scalar operations
    // in O(1), without allocation, cloning, branching, or parameter lookup.
    let dist_term = distance - r0;
    let dist_term2 = dist_term * dist_term;
    let c1 = 143.9325;
    let cs = -2.0;
    let c3 = 7.0 / 12.0;

    0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2)
}
pub(super) fn calc_angle_bend_energy(theta0: f64, ka: f64, is_linear: bool, cos_theta: f64) -> f64 {
    const M_PI: f64 = std::f64::consts::PI;
    const DEG2RAD: f64 = M_PI / 180.0;
    const RAD2DEG: f64 = 180.0 / M_PI;
    const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/AngleBend.cpp:43-57:
    // RDKit❗✔️: double calcAngleBendEnergy(const double theta0, const double ka, bool isLinear,
    // RDKit❗✔️:                            const double cosTheta) {
    // RDKit❗✔️:   double angle = RAD2DEG * acos(cosTheta) - theta0;
    // RDKit❗✔️:   double const cb = -0.006981317;
    // RDKit❗✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    // RDKit❗✔️:   double res = 0.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   if (isLinear) {
    // RDKit❗✔️:     res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cosTheta);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️: #ifndef M_PI
    // RDKit❗✔️: #define M_PI 3.14159265358979323846
    // RDKit❗✔️: #endif
    // RDKit❗✔️: constexpr double DEG2RAD = M_PI / 180.0;
    // RDKit❗✔️: constexpr double RAD2DEG = 180.0 / M_PI;
    // RDKit❗✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
    // Behavior review: compute angle via acos before the linear branch, do not
    // clamp cosine, and preserve both source formulas and operation grouping.
    // Complexity review: one libm acos and fixed scalar arithmetic per call;
    // O(1), without allocation, cloning, searching, or parameter lookup.
    let angle = RAD2DEG * cos_theta.acos() - theta0;
    let cb = -0.006981317;
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    let mut res = 0.0;

    if is_linear {
        res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cos_theta);
    } else {
        res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
    }

    res
}

#[cfg(test)]
mod tests {
    use super::{calc_angle_bend_energy, calc_bond_stretch_energy};
    const THETA0_BITS: [u64; 3] = [
        0x0000_0000_0000_0000,
        0x405b_5e14_7ae1_47ae,
        0x4066_8000_0000_0000,
    ];
    const KA_BITS: [u64; 3] = [
        0x8000_0000_0000_0000,
        0x3fe0_0000_0000_0000,
        0x4010_0000_0000_0000,
    ];
    const COS_THETA_BITS: [u64; 5] = [
        0xbff0_0000_0000_0000,
        0xbfe0_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x3fe0_0000_0000_0000,
        0x3ff0_0000_0000_0000,
    ];
    // Literal input and result bits copied from the frozen pinned C++ source
    // driver. The test never invokes or reads the native oracle.
    const SOURCE_ANGLE_BITS: [[u64; 5]; 94] = [
        [
            0x0000000000000000,
            0x8000000000000000,
            0,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            0,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            0,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            0,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            0,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            1,
            0xbff0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            1,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            1,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            1,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x8000000000000000,
            1,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            0,
            0xbff0000000000000,
            0xc056c9149a24c3de,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            0,
            0xbfe0000000000000,
            0x40399bb3e84da6bc,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            0,
            0x0000000000000000,
            0x40407ff50c89f34f,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            0,
            0x3fe0000000000000,
            0x4036ee54e3539c33,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            0,
            0x3ff0000000000000,
            0x0000000000000000,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            1,
            0xbfe0000000000000,
            0x4041fdd70a3d70a4,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            1,
            0x0000000000000000,
            0x4051fdd70a3d70a4,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            1,
            0x3fe0000000000000,
            0x405afcc28f5c28f6,
        ],
        [
            0x0000000000000000,
            0x3fe0000000000000,
            1,
            0x3ff0000000000000,
            0x4061fdd70a3d70a4,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            0,
            0xbff0000000000000,
            0xc086c9149a24c3de,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            0,
            0xbfe0000000000000,
            0x40699bb3e84da6bc,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            0,
            0x0000000000000000,
            0x40707ff50c89f34f,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            0,
            0x3fe0000000000000,
            0x4066ee54e3539c33,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            0,
            0x3ff0000000000000,
            0x0000000000000000,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            1,
            0xbfe0000000000000,
            0x4071fdd70a3d70a4,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            1,
            0x0000000000000000,
            0x4081fdd70a3d70a4,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            1,
            0x3fe0000000000000,
            0x408afcc28f5c28f6,
        ],
        [
            0x0000000000000000,
            0x4010000000000000,
            1,
            0x3ff0000000000000,
            0x4091fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            0,
            0xbff0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            0,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            0,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            0,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            0,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            1,
            0xbff0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            1,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            1,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            1,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x8000000000000000,
            1,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0xbff0000000000000,
            0x403bad7c0d86db88,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0xbfe0000000000000,
            0x3ff20436f0c4fc6a,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0x0000000000000000,
            0x4012e135968b18ea,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0x3fe0000000000000,
            0x40420b6c64adf65f,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0x3ff0000000000000,
            0x406cf7b5727edc32,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0xbfe0000000000000,
            0x4041fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0x0000000000000000,
            0x4051fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0x3fe0000000000000,
            0x405afcc28f5c28f6,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0x3ff0000000000000,
            0x4061fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            0,
            0xbff0000000000000,
            0x406bad7c0d86db88,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            0,
            0xbfe0000000000000,
            0x40220436f0c4fc6a,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            0,
            0x0000000000000000,
            0x4042e135968b18ea,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            0,
            0x3fe0000000000000,
            0x40720b6c64adf65f,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            0,
            0x3ff0000000000000,
            0x409cf7b5727edc32,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            1,
            0xbfe0000000000000,
            0x4071fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            1,
            0x0000000000000000,
            0x4081fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            1,
            0x3fe0000000000000,
            0x408afcc28f5c28f6,
        ],
        [
            0x405b5e147ae147ae,
            0x4010000000000000,
            1,
            0x3ff0000000000000,
            0x4091fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            0,
            0xbff0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            0,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            0,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            0,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            0,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            1,
            0xbff0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            1,
            0xbfe0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            1,
            0x0000000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            1,
            0x3fe0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x8000000000000000,
            1,
            0x3ff0000000000000,
            0x8000000000000000,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            0,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            0,
            0xbfe0000000000000,
            0x404bfe925aea0098,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            0,
            0x0000000000000000,
            0x4062123ceff0a772,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            0,
            0x3fe0000000000000,
            0x4072212327c50cef,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            0,
            0x3ff0000000000000,
            0x40890b5cc657bcc2,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            1,
            0xbfe0000000000000,
            0x4041fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            1,
            0x0000000000000000,
            0x4051fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            1,
            0x3fe0000000000000,
            0x405afcc28f5c28f6,
        ],
        [
            0x4066800000000000,
            0x3fe0000000000000,
            1,
            0x3ff0000000000000,
            0x4061fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            0,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            0,
            0xbfe0000000000000,
            0x407bfe925aea0098,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            0,
            0x0000000000000000,
            0x4092123ceff0a772,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            0,
            0x3fe0000000000000,
            0x40a2212327c50cef,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            0,
            0x3ff0000000000000,
            0x40b90b5cc657bcc2,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            1,
            0xbff0000000000000,
            0x0000000000000000,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            1,
            0xbfe0000000000000,
            0x4071fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            1,
            0x0000000000000000,
            0x4081fdd70a3d70a4,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            1,
            0x3fe0000000000000,
            0x408afcc28f5c28f6,
        ],
        [
            0x4066800000000000,
            0x4010000000000000,
            1,
            0x3ff0000000000000,
            0x4091fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0xbff0000000000001,
            0x7ff8000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0xbff0000000000001,
            0xbd11fdd70a3d70a4,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            0,
            0x3ff0000000000001,
            0x7ff8000000000000,
        ],
        [
            0x405b5e147ae147ae,
            0x3fe0000000000000,
            1,
            0x3ff0000000000001,
            0x4061fdd70a3d70a4,
        ],
    ];

    const R0_BITS: [u64; 3] = [
        0x0000_0000_0000_0000,
        0x3ff4_0000_0000_0000,
        0x4000_0000_0000_0000,
    ];
    const KB_BITS: [u64; 3] = [
        0x8000_0000_0000_0000,
        0x3fe0_0000_0000_0000,
        0x4010_0000_0000_0000,
    ];
    const DISTANCE_BITS: [u64; 5] = [
        0xbff0_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x3ff4_0000_0000_0000,
        0x4000_0000_0000_0000,
        0x4014_0000_0000_0000,
    ];

    // Input bits and source-expression output bits from the independent
    // pinned C++20 source-only driver; rows are in the frozen Cartesian order.
    const SOURCE_BITS: [[u64; 4]; 45] = [
        [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x4067_fd1e_b851_eb86,
        ],
        [
            0x0000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x405e_2961_0000_0001,
        ],
        [
            0x0000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x408c_7c94_7ae1_47af,
        ],
        [
            0x0000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x40e5_ab66_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x4097_fd1e_b851_eb86,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x408e_2961_0000_0001,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x40bc_7c94_7ae1_47af,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x4115_ab66_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x40a8_a372_c051_eb86,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x4079_1c3c_4000_0001,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x4030_7206_8f5c_28f6,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x40ca_013a_c200_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x4010_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x40d8_a372_c051_eb86,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x40a9_1c3c_4000_0001,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x4060_7206_8f5c_28f6,
        ],
        [
            0x3ff4_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x40fa_013a_c200_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x40c1_b5df_ae14_7ae1,
        ],
        [
            0x4000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x40a0_1e10_a3d7_0a3e,
        ],
        [
            0x4000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x4053_4aaf_147a_e147,
        ],
        [
            0x4000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x40b4_3d91_eb85_1eb8,
        ],
        [
            0x4000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x40f1_b5df_ae14_7ae1,
        ],
        [
            0x4000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x40d0_1e10_a3d7_0a3e,
        ],
        [
            0x4000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x4083_4aaf_147a_e147,
        ],
        [
            0x4000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x4000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4014_0000_0000_0000,
            0x40e4_3d91_eb85_1eb8,
        ],
    ];

    #[test]
    fn mmff_bonded_bond_90_calls_match_source_bits_and_preserve_inputs() {
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;
        let mut base_calls = 0;

        for (row_index, row) in SOURCE_BITS.iter().enumerate() {
            let expected_inputs = [
                R0_BITS[row_index / 15],
                KB_BITS[(row_index % 15) / 5],
                DISTANCE_BITS[row_index % 5],
            ];
            if [row[0], row[1], row[2]] != expected_inputs {
                discrepancies.push(format!(
                    "source row {row_index}: frozen input bits differ from fixed axes"
                ));
            }

            for repeat in 0..2 {
                let frozen_row = *row;
                let row_address_before = row as *const [u64; 4] as usize;
                let case = format!("source row {row_index} repeat {repeat}");
                let r0 = f64::from_bits(frozen_row[0]);
                let kb = f64::from_bits(frozen_row[1]);
                let distance = f64::from_bits(frozen_row[2]);
                let input_before = [r0.to_bits(), kb.to_bits(), distance.to_bits()];
                if input_before != expected_inputs {
                    discrepancies.push(format!(
                        "{case}: decoded scalar inputs differ from fixed axes"
                    ));
                }

                let actual = calc_bond_stretch_energy(r0, kb, distance);
                let actual_bits = actual.to_bits();
                actual_calls += 1;
                base_calls += 1;
                if repeat == 0 {
                    base_cells += 1;
                }

                let input_after = [r0.to_bits(), kb.to_bits(), distance.to_bits()];
                let row_after = *row;
                let row_address_after = row as *const [u64; 4] as usize;
                if input_after != input_before {
                    discrepancies.push(format!("{case}: scalar input bits changed"));
                }
                if row_after != frozen_row || row_address_after != row_address_before {
                    discrepancies.push(format!(
                        "{case}: immutable source row value or address changed"
                    ));
                }
                if actual_bits != frozen_row[3] {
                    discrepancies.push(format!(
                        "{case}: output bits {actual_bits:016x}, expected {:016x}",
                        frozen_row[3]
                    ));
                }
            }
        }

        if actual_calls != 90 {
            discrepancies.push(format!("expected 90 actual bond calls, got {actual_calls}"));
        }
        if base_cells != 45 {
            discrepancies.push(format!("expected 45 repeated bond cells, got {base_cells}"));
        }
        if base_calls != 90 {
            discrepancies.push(format!("expected 90 repeated bond calls, got {base_calls}"));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF bond stretch energy discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_bonded_angle_184_calls_match_source_bits_and_preserve_inputs() {
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;
        let mut base_calls = 0;
        let mut supplement_cells = 0;
        let mut supplement_calls = 0;

        if std::f64::consts::PI.to_bits() != 0x4009_21fb_5444_2d18 {
            discrepancies.push(String::from(
                "Rust PI bits differ from the frozen source M_PI literal",
            ));
        }

        for (row_index, row) in SOURCE_ANGLE_BITS.iter().enumerate() {
            let expected_inputs = if row_index < 90 {
                [
                    THETA0_BITS[row_index / 30],
                    KA_BITS[(row_index % 30) / 10],
                    ((row_index % 10) / 5) as u64,
                    COS_THETA_BITS[row_index % 5],
                ]
            } else {
                let supplement_index = row_index - 90;
                [
                    THETA0_BITS[1],
                    KA_BITS[1],
                    (supplement_index % 2) as u64,
                    if supplement_index < 2 {
                        0xbff0_0000_0000_0001
                    } else {
                        0x3ff0_0000_0000_0001
                    },
                ]
            };
            if [row[0], row[1], row[2], row[3]] != expected_inputs {
                discrepancies.push(format!(
                    "source angle row {row_index}: frozen input bits differ from fixed axes"
                ));
            }

            let repeats = if row_index < 90 { 2 } else { 1 };
            for repeat in 0..repeats {
                let frozen_row = *row;
                let row_address_before = row as *const [u64; 5] as usize;
                let case = format!("source angle row {row_index} repeat {repeat}");
                let theta0 = f64::from_bits(frozen_row[0]);
                let ka = f64::from_bits(frozen_row[1]);
                let is_linear = match frozen_row[2] {
                    0 => false,
                    1 => true,
                    value => {
                        discrepancies
                            .push(format!("{case}: invalid frozen is_linear value {value}"));
                        value != 0
                    }
                };
                let cos_theta = f64::from_bits(frozen_row[3]);
                let input_before = [
                    theta0.to_bits(),
                    ka.to_bits(),
                    is_linear as u64,
                    cos_theta.to_bits(),
                ];
                if input_before != expected_inputs {
                    discrepancies.push(format!(
                        "{case}: decoded scalar inputs differ from fixed axes"
                    ));
                }

                let actual = calc_angle_bend_energy(theta0, ka, is_linear, cos_theta);
                let actual_bits = actual.to_bits();
                actual_calls += 1;
                if row_index < 90 {
                    base_calls += 1;
                    if repeat == 0 {
                        base_cells += 1;
                    }
                } else {
                    supplement_calls += 1;
                    if repeat == 0 {
                        supplement_cells += 1;
                    }
                }

                let input_after = [
                    theta0.to_bits(),
                    ka.to_bits(),
                    is_linear as u64,
                    cos_theta.to_bits(),
                ];
                let row_after = *row;
                let row_address_after = row as *const [u64; 5] as usize;
                if input_after != input_before {
                    discrepancies.push(format!("{case}: scalar input bits changed"));
                }
                if row_after != frozen_row || row_address_after != row_address_before {
                    discrepancies.push(format!(
                        "{case}: immutable source row value or address changed"
                    ));
                }
                if actual_bits != frozen_row[4] {
                    discrepancies.push(format!(
                        "{case}: output bits {actual_bits:016x}, expected {:016x}",
                        frozen_row[4]
                    ));
                }
            }
        }

        if base_cells != 90 {
            discrepancies.push(format!("expected 90 base angle cells, got {base_cells}"));
        }
        if base_calls != 180 {
            discrepancies.push(format!(
                "expected 180 repeated base angle calls, got {base_calls}"
            ));
        }
        if supplement_cells != 4 || supplement_calls != 4 {
            discrepancies.push(format!(
                "expected 4 once-only supplement cells/calls, got {supplement_cells}/{supplement_calls}"
            ));
        }
        if actual_calls != 184 {
            discrepancies.push(format!(
                "expected 184 actual angle calls, got {actual_calls}"
            ));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF angle bend energy discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }
}
