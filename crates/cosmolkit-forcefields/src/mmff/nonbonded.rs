use super::params::{MmffVdw, MmffVdwCollection};

pub(super) fn calc_unscaled_vdw_minimum(
    collection: &MmffVdwCollection,
    i_atom: &MmffVdw,
    j_atom: &MmffVdw,
) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Nonbonded.cpp:20-32:
    // RDKit✔️✔️: double calcUnscaledVdWMinimum(const MMFFVdWCollection *mmffVdW,
    // RDKit✔️✔️:                               const MMFFVdW *mmffVdWParamsIAtom,
    // RDKit✔️✔️:                               const MMFFVdW *mmffVdWParamsJAtom) {
    // RDKit✔️✔️:   double gamma_ij = (mmffVdWParamsIAtom->R_star - mmffVdWParamsJAtom->R_star) /
    // RDKit✔️✔️:                     (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return (0.5 * (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star) *
    // RDKit✔️✔️:           (1.0 +
    // RDKit✔️✔️:            (((mmffVdWParamsIAtom->DA == 'D') || (mmffVdWParamsJAtom->DA == 'D'))
    // RDKit✔️✔️:                 ? 0.0
    // RDKit✔️✔️:                 : mmffVdW->B *
    // RDKit✔️✔️:                       (1.0 - exp(-(mmffVdW->Beta) * gamma_ij * gamma_ij)))));
    // RDKit✔️✔️: }
    // Behavior review: retain the source gamma, donor short-circuit, unary
    // minus, exponential and multiplication grouping for these borrowed rows.
    // Complexity review: direct field reads and one conditional exponential
    // are O(1), with no allocation, clone, scan, or collection lookup.
    let gamma_ij = (i_atom.r_star - j_atom.r_star) / (i_atom.r_star + j_atom.r_star);

    0.5 * (i_atom.r_star + j_atom.r_star)
        * (1.0
            + if i_atom.da == b'D' || j_atom.da == b'D' {
                0.0
            } else {
                collection.b * (1.0 - (-(collection.beta) * gamma_ij * gamma_ij).exp())
            })
}

pub(super) fn calc_unscaled_vdw_well_depth(
    r_star_ij: f64,
    i_atom: &MmffVdw,
    j_atom: &MmffVdw,
) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Nonbonded.cpp:34-45:
    // RDKit✔️✔️: double calcUnscaledVdWWellDepth(double R_star_ij,
    // RDKit✔️✔️:                                 const MMFFVdW *mmffVdWParamsIAtom,
    // RDKit✔️✔️:                                 const MMFFVdW *mmffVdWParamsJAtom) {
    // RDKit✔️✔️:   double R_star_ij2 = R_star_ij * R_star_ij;
    // RDKit✔️✔️:   double const c4 = 181.16;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return (c4 * mmffVdWParamsIAtom->G_i * mmffVdWParamsJAtom->G_i *
    // RDKit✔️✔️:           mmffVdWParamsIAtom->alpha_i * mmffVdWParamsJAtom->alpha_i /
    // RDKit✔️✔️:           ((sqrt(mmffVdWParamsIAtom->alpha_i / mmffVdWParamsIAtom->N_i) +
    // RDKit✔️✔️:             sqrt(mmffVdWParamsJAtom->alpha_i / mmffVdWParamsJAtom->N_i)) *
    // RDKit✔️✔️:            R_star_ij2 * R_star_ij2 * R_star_ij2));
    // RDKit✔️✔️: }
    // Behavior review: square the independent supplied radius once; retain
    // the left-associated numerator and source denominator multiplication.
    // Complexity review: fixed field reads and arithmetic/square roots are
    // O(1), with no allocation, clone, scan, or collection lookup.
    let r_star_ij2 = r_star_ij * r_star_ij;
    let c4 = 181.16;

    (c4 * i_atom.g_i * j_atom.g_i * i_atom.alpha_i * j_atom.alpha_i)
        / (((i_atom.alpha_i / i_atom.n_i).sqrt() + (j_atom.alpha_i / j_atom.n_i).sqrt())
            * r_star_ij2
            * r_star_ij2
            * r_star_ij2)
}

pub(super) fn scale_vdw_params(
    r_star_ij: &mut f64,
    well_depth: &mut f64,
    collection: &MmffVdwCollection,
    i_atom: &MmffVdw,
    j_atom: &MmffVdw,
) {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Nonbonded.cpp:66-75:
    // RDKit✔️✔️: void scaleVdWParams(double &R_star_ij, double &wellDepth,
    // RDKit✔️✔️:                     const MMFFVdWCollection *mmffVdW,
    // RDKit✔️✔️:                     const MMFFVdW *mmffVdWParamsIAtom,
    // RDKit✔️✔️:                     const MMFFVdW *mmffVdWParamsJAtom) {
    // RDKit✔️✔️:   if (((mmffVdWParamsIAtom->DA == 'D') && (mmffVdWParamsJAtom->DA == 'A')) ||
    // RDKit✔️✔️:       ((mmffVdWParamsIAtom->DA == 'A') && (mmffVdWParamsJAtom->DA == 'D'))) {
    // RDKit✔️✔️:     R_star_ij *= mmffVdW->DARAD;
    // RDKit✔️✔️:     wellDepth *= mmffVdW->DAEPS;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: mutate only ordered donor/acceptor pairs, radius first.
    // Complexity review: two DA comparisons and at most two scalar multiplies;
    // fixed O(1) work, no allocation, clone, scan, or collection lookup.
    if i_atom.da == b'D' && j_atom.da == b'A' || i_atom.da == b'A' && j_atom.da == b'D' {
        *r_star_ij *= collection.darad;
        *well_depth *= collection.daeps;
    }
}

pub(super) fn calc_vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Nonbonded.cpp:47-64:
    // RDKit✔️✔️: double calcVdWEnergy(const double dist, const double R_star_ij,
    // RDKit✔️✔️:                      const double wellDepth) {
    // RDKit✔️✔️:   double const vdw1 = 1.07;
    // RDKit✔️✔️:   double const vdw1m1 = vdw1 - 1.0;
    // RDKit✔️✔️:   double const vdw2 = 1.12;
    // RDKit✔️✔️:   double const vdw2m1 = vdw2 - 1.0;
    // RDKit✔️✔️:   double dist2 = dist * dist;
    // RDKit✔️✔️:   double dist7 = dist2 * dist2 * dist2 * dist;
    // RDKit✔️✔️:   double aTerm = vdw1 * R_star_ij / (dist + vdw1m1 * R_star_ij);
    // RDKit✔️✔️:   double aTerm2 = aTerm * aTerm;
    // RDKit✔️✔️:   double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
    // RDKit✔️✔️:   double R_star_ij2 = R_star_ij * R_star_ij;
    // RDKit✔️✔️:   double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
    // RDKit✔️✔️:   double bTerm = vdw2 * R_star_ij7 / (dist7 + vdw2m1 * R_star_ij7) - 2.0;
    // RDKit✔️✔️:   double res = wellDepth * aTerm7 * bTerm;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Behavior review: preserve the complete source expression and every
    // intermediate/product grouping, including IEEE behavior for signed zero
    // and the frozen infinite-distance inputs.
    // Complexity review: fixed scalar arithmetic only; O(1) work with no
    // allocation, clone, table access, scan, or branch.
    let vdw1 = 1.07;
    let vdw1m1 = vdw1 - 1.0;
    let vdw2 = 1.12;
    let vdw2m1 = vdw2 - 1.0;
    let dist2 = dist * dist;
    let dist7 = dist2 * dist2 * dist2 * dist;
    let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
    let a_term2 = a_term * a_term;
    let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
    let r_star_ij2 = r_star_ij * r_star_ij;
    let r_star_ij7 = r_star_ij2 * r_star_ij2 * r_star_ij2 * r_star_ij;
    let b_term = vdw2 * r_star_ij7 / (dist7 + vdw2m1 * r_star_ij7) - 2.0;
    let res = well_depth * a_term7 * b_term;

    res
}

pub(super) fn calc_ele_energy(
    _idx1: u32,
    _idx2: u32,
    dist: f64,
    charge_term: f64,
    diel_model: u8,
    is_1_4: bool,
) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h:67-70:
    // RDKit✔️✔️: enum {
    // RDKit✔️✔️:   CONSTANT = 1,
    // RDKit✔️✔️:   DISTANCE = 2
    // RDKit✔️✔️: };
    const DISTANCE: u8 = 2;

    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Nonbonded.cpp:77-86:
    // RDKit✔️✔️: double calcEleEnergy(unsigned int, unsigned int, double dist, double chargeTerm,
    // RDKit✔️✔️:                      std::uint8_t dielModel, bool is1_4) {
    // RDKit✔️✔️:   double corr_dist = dist + 0.05;
    // RDKit✔️✔️:   double const diel = 332.0716;
    // RDKit✔️✔️:   double const sc1_4 = 0.75;
    // RDKit✔️✔️:   if (dielModel == RDKit::MMFF::DISTANCE) {
    // RDKit✔️✔️:     corr_dist *= corr_dist;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return (diel * chargeTerm / corr_dist * (is1_4 ? sc1_4 : 1.0));
    // RDKit✔️✔️: }
    // Behavior review: retain the ignored indices, exact byte comparison,
    // source corrected-distance square, and left-associated Coulomb expression
    // followed by the conditional 1-4 factor.
    // Complexity review: fixed scalar arithmetic and one source branch are
    // O(1), with no allocation, clone, lookup, scan, or fallback.
    let mut corr_dist = dist + 0.05;
    let diel = 332.0716;
    let sc1_4 = 0.75;
    if diel_model == DISTANCE {
        corr_dist *= corr_dist;
    }
    (diel * charge_term / corr_dist) * (if is_1_4 { sc1_4 } else { 1.0 })
}

#[cfg(test)]
mod tests {
    use super::{calc_ele_energy, calc_vdw_energy};

    const BASE_DISTANCE_BITS: [u64; 3] =
        [0x0000000000000000, 0x3fe8000000000000, 0x4020000000000000];
    const BASE_CHARGE_BITS: [u64; 3] = [0x8000000000000000, 0x3fe4000000000000, 0xbff8000000000000];
    const MODEL_BYTES: [u8; 4] = [1, 2, 0, 255];
    const BOUNDARY_DISTANCE_BITS: [u64; 2] = [0xbf9999999999999a, 0x7ff0000000000000];
    const BOUNDARY_MODEL_BYTES: [u8; 2] = [1, 2];

    // idx1, idx2, distance bits, charge bits, dielectric byte, 1-4 flag,
    // and independent source-expression output bits.
    const ELE_SOURCE_BITS: [[u64; 7]; 88] = [
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x01,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x01,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x02,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x02,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x00,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0x00,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0xff,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x8000000000000000,
            0xff,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x01,
            0,
            0x40b036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x01,
            1,
            0x40a85257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x02,
            0,
            0x40f4449e66666665,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x02,
            1,
            0x40ee66ed99999998,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x00,
            0,
            0x40b036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0x00,
            1,
            0x40a85257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0xff,
            0,
            0x40b036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0x3fe4000000000000,
            0xff,
            1,
            0x40a85257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x01,
            0,
            0xc0c37512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x01,
            1,
            0xc0bd2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x02,
            0,
            0xc1085257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x02,
            1,
            0xc1023dc1c28f5c28,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x00,
            0,
            0xc0c37512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0x00,
            1,
            0xc0bd2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0xff,
            0,
            0xc0c37512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x0000000000000000,
            0xbff8000000000000,
            0xff,
            1,
            0xc0bd2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x01,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x01,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x02,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x02,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x00,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0x00,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0xff,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x8000000000000000,
            0xff,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x01,
            0,
            0x407036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x01,
            1,
            0x40685257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x02,
            0,
            0x4074449e66666665,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x02,
            1,
            0x406e66ed99999998,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x00,
            0,
            0x407036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0x00,
            1,
            0x40685257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0xff,
            0,
            0x407036e51eb851eb,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0x3fe4000000000000,
            0xff,
            1,
            0x40685257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x01,
            0,
            0xc0837512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x01,
            1,
            0xc07d2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x02,
            0,
            0xc0885257ae147ae0,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x02,
            1,
            0xc0823dc1c28f5c28,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x00,
            0,
            0xc0837512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0x00,
            1,
            0xc07d2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0xff,
            0,
            0xc0837512f1a9fbe7,
        ],
        [
            0x00000000,
            0x00000001,
            0x3fe8000000000000,
            0xbff8000000000000,
            0xff,
            1,
            0xc07d2f9c6a7ef9da,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x01,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x01,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x02,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x02,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x00,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0x00,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0xff,
            0,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x8000000000000000,
            0xff,
            1,
            0x8000000000000000,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x01,
            0,
            0x4039c82e4d77c372,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x01,
            1,
            0x40335622ba19d296,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x02,
            0,
            0x40099f2f9ae652ed,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x02,
            1,
            0x40033763b42cbe32,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x00,
            0,
            0x4039c82e4d77c372,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0x00,
            1,
            0x40335622ba19d296,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0xff,
            0,
            0x4039c82e4d77c372,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0x3fe4000000000000,
            0xff,
            1,
            0x40335622ba19d296,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x01,
            0,
            0xc04ef037902950f0,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x01,
            1,
            0xc0473429ac1efcb4,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x02,
            0,
            0xc01ebf05ed146383,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x02,
            1,
            0xc0170f4471cf4aa2,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x00,
            0,
            0xc04ef037902950f0,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0x00,
            1,
            0xc0473429ac1efcb4,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0xff,
            0,
            0xc04ef037902950f0,
        ],
        [
            0x00000000,
            0x00000001,
            0x4020000000000000,
            0xbff8000000000000,
            0xff,
            1,
            0xc0473429ac1efcb4,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x01,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x01,
            1,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x02,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x02,
            1,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x00,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0x00,
            1,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0xff,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x3fe8000000000000,
            0x0000000000000000,
            0xff,
            1,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0xbf9999999999999a,
            0x3fe4000000000000,
            0x01,
            0,
            0x40c036e51eb851eb,
        ],
        [
            0x00000011,
            0x0000001d,
            0xbf9999999999999a,
            0x3fe4000000000000,
            0x01,
            1,
            0x40b85257ae147ae0,
        ],
        [
            0x00000011,
            0x0000001d,
            0xbf9999999999999a,
            0x3fe4000000000000,
            0x02,
            0,
            0x4114449e66666665,
        ],
        [
            0x00000011,
            0x0000001d,
            0xbf9999999999999a,
            0x3fe4000000000000,
            0x02,
            1,
            0x410e66ed99999998,
        ],
        [
            0x00000011,
            0x0000001d,
            0x7ff0000000000000,
            0x3fe4000000000000,
            0x01,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x7ff0000000000000,
            0x3fe4000000000000,
            0x01,
            1,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x7ff0000000000000,
            0x3fe4000000000000,
            0x02,
            0,
            0x0000000000000000,
        ],
        [
            0x00000011,
            0x0000001d,
            0x7ff0000000000000,
            0x3fe4000000000000,
            0x02,
            1,
            0x0000000000000000,
        ],
    ];

    #[test]
    fn mmff_nonbonded_energy_ele_160_calls_match_source_bits_and_preserve_inputs() {
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;
        let mut base_calls = 0;
        let mut supplementary_rows = 0;
        let mut supplementary_calls = 0;

        for (row_index, row) in ELE_SOURCE_BITS.iter().enumerate() {
            let repetitions = if row_index < 72 { 2 } else { 1 };

            for repeat in 0..repetitions {
                let frozen_row = *row;
                let case = format!("source row {row_index} repeat {repeat}");
                let row_address_before = row as *const [u64; 7] as usize;
                let expected_frozen_inputs = if row_index < 72 {
                    [
                        0,
                        1,
                        BASE_DISTANCE_BITS[row_index / 24],
                        BASE_CHARGE_BITS[(row_index % 24) / 8],
                        u64::from(MODEL_BYTES[(row_index % 8) / 2]),
                        (row_index % 2) as u64,
                    ]
                } else if row_index < 80 {
                    let local_index = row_index - 72;
                    [
                        17,
                        29,
                        0x3fe8000000000000,
                        0,
                        u64::from(MODEL_BYTES[local_index / 2]),
                        (local_index % 2) as u64,
                    ]
                } else {
                    let local_index = row_index - 80;
                    [
                        17,
                        29,
                        BOUNDARY_DISTANCE_BITS[local_index / 4],
                        0x3fe4000000000000,
                        u64::from(BOUNDARY_MODEL_BYTES[(local_index % 4) / 2]),
                        (local_index % 2) as u64,
                    ]
                };
                let frozen_inputs = [
                    frozen_row[0],
                    frozen_row[1],
                    frozen_row[2],
                    frozen_row[3],
                    frozen_row[4],
                    frozen_row[5],
                ];
                if frozen_inputs != expected_frozen_inputs {
                    discrepancies.push(format!(
                        "{case}: literal input row differs from fixed matrix"
                    ));
                }

                let expected_indexes = if row_index < 72 && repeat == 1 {
                    [u32::MAX, u32::MAX]
                } else {
                    [frozen_row[0] as u32, frozen_row[1] as u32]
                };
                let idx1 = expected_indexes[0];
                let idx2 = expected_indexes[1];
                let dist = f64::from_bits(frozen_row[2]);
                let charge_term = f64::from_bits(frozen_row[3]);
                let diel_model = frozen_row[4] as u8;
                let is_1_4 = frozen_row[5] != 0;
                let input_before = [
                    u64::from(idx1),
                    u64::from(idx2),
                    dist.to_bits(),
                    charge_term.to_bits(),
                    u64::from(diel_model),
                    is_1_4 as u64,
                ];
                let expected_input = [
                    u64::from(expected_indexes[0]),
                    u64::from(expected_indexes[1]),
                    frozen_row[2],
                    frozen_row[3],
                    frozen_row[4],
                    frozen_row[5],
                ];
                if input_before != expected_input {
                    discrepancies.push(format!("{case}: decoded input differs from frozen row"));
                }

                let actual = calc_ele_energy(idx1, idx2, dist, charge_term, diel_model, is_1_4);
                let actual_bits = actual.to_bits();
                actual_calls += 1;
                if row_index < 72 {
                    base_calls += 1;
                    if repeat == 0 {
                        base_cells += 1;
                    }
                } else {
                    supplementary_calls += 1;
                    if repeat == 0 {
                        supplementary_rows += 1;
                    }
                }

                let input_after = [
                    u64::from(idx1),
                    u64::from(idx2),
                    dist.to_bits(),
                    charge_term.to_bits(),
                    u64::from(diel_model),
                    is_1_4 as u64,
                ];
                let row_after = *row;
                let row_address_after = row as *const [u64; 7] as usize;
                if input_after != input_before {
                    discrepancies.push(format!(
                        "{case}: scalar input bits, indexes, or flag changed"
                    ));
                }
                if row_after != frozen_row || row_address_after != row_address_before {
                    discrepancies.push(format!(
                        "{case}: immutable source row value or address changed"
                    ));
                }
                if actual_bits != frozen_row[6] {
                    discrepancies.push(format!(
                        "{case}: output bits {actual_bits:016x}, expected {:016x}",
                        frozen_row[6]
                    ));
                }
            }
        }

        if actual_calls != 160 {
            discrepancies.push(format!(
                "expected 160 actual electrostatic calls, got {actual_calls}"
            ));
        }
        if base_cells != 72 {
            discrepancies.push(format!(
                "expected 72 base electrostatic cells, got {base_cells}"
            ));
        }
        if base_calls != 144 {
            discrepancies.push(format!(
                "expected 144 repeated base calls, got {base_calls}"
            ));
        }
        if supplementary_rows != 16 {
            discrepancies.push(format!(
                "expected 16 supplementary rows, got {supplementary_rows}"
            ));
        }
        if supplementary_calls != 16 {
            discrepancies.push(format!(
                "expected 16 supplementary calls, got {supplementary_calls}"
            ));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF electrostatic energy discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    static SOURCE_BITS: [[u64; 4]; 24] = [
        [
            0x0000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x41c5_4e95_9cb2_1fb8,
        ],
        [
            0x0000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x41c5_4e95_9cb2_1fb3,
        ],
        [
            0x0000_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3fc0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4135_e4fd_f510_0165,
        ],
        [
            0x3fc0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3fc0_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4189_bcbe_bfe3_6e75,
        ],
        [
            0x3fc0_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4007_c877_770a_01e3,
        ],
        [
            0x3ff0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x3ff0_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x40d0_b936_d0b6_b90c,
        ],
        [
            0x3ff0_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x400c_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbf50_6585_3174_fc36,
        ],
        [
            0x400c_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x400c_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3fcc_2d9b_ffdc_b6fd,
        ],
        [
            0x400c_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x4028_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbe8b_4252_2b1a_160c,
        ],
        [
            0x4028_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x4028_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbf44_7794_9046_6e15,
        ],
        [
            0x4028_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x8000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x41c5_4e95_9cb2_1fb8,
        ],
        [
            0x8000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x7ff0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x8000_0000_0000_0000,
        ],
        [
            0x7ff0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
    ];

    #[test]
    fn mmff_nonbonded_energy_vdw_44_calls_match_source_bits_and_preserve_inputs() {
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_calls = 0;
        let mut supplementary_calls = 0;

        for (row_index, row) in SOURCE_BITS.iter().enumerate() {
            let repetitions = if row_index < 20 { 2 } else { 1 };
            for repeat in 0..repetitions {
                let case = format!("source row {row_index} repeat {repeat}");
                let row_before = *row;
                let row_address_before = row as *const [u64; 4] as usize;
                let dist = f64::from_bits(row_before[0]);
                let r_star_ij = f64::from_bits(row_before[1]);
                let well_depth = f64::from_bits(row_before[2]);
                let input_bits_before = [dist.to_bits(), r_star_ij.to_bits(), well_depth.to_bits()];
                if input_bits_before != [row_before[0], row_before[1], row_before[2]] {
                    discrepancies.push(format!("{case}: decoded inputs differ from frozen row"));
                }

                // The source function accepts scalar values and returns only
                // the energy; no collection, molecule, or cache state exists.
                let actual = calc_vdw_energy(dist, r_star_ij, well_depth);
                let input_bits_after = [dist.to_bits(), r_star_ij.to_bits(), well_depth.to_bits()];
                let row_after = *row;
                let row_address_after = row as *const [u64; 4] as usize;
                let actual_bits = actual.to_bits();
                actual_calls += 1;

                if input_bits_after != input_bits_before
                    || row_after != row_before
                    || row_address_after != row_address_before
                {
                    discrepancies.push(format!("{case}: scalar input or literal row changed"));
                }
                if actual_bits != row_before[3] {
                    discrepancies.push(format!(
                        "{case}: output bits {actual_bits:016x}, expected {:016x}",
                        row_before[3]
                    ));
                }
                if row_index < 20 {
                    base_calls += 1;
                } else {
                    supplementary_calls += 1;
                }
            }
        }

        if actual_calls != 44 {
            discrepancies.push(format!("expected 44 actual VdW calls, got {actual_calls}"));
        }
        if base_calls != 40 {
            discrepancies.push(format!("expected 40 repeated base calls, got {base_calls}"));
        }
        if supplementary_calls != 4 {
            discrepancies.push(format!(
                "expected 4 signed-zero/infinite-distance controls, got {supplementary_calls}"
            ));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW energy discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }
}
