//! BIO-owned crystallographic coordinate primitives.

use super::{BioCrystalInfo, multiply_matrix_vector};
use crate::BioAsu;

/// A nearest-image result with Gemmi's squared distance, periodic shift and
/// one-based symmetry index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BioNearestImage {
    dist_sq: f64,
    pbc_shift: [i32; 3],
    sym_idx: i32,
}

impl BioNearestImage {
    pub(super) const fn from_parts(dist_sq: f64, pbc_shift: [i32; 3], sym_idx: i32) -> Self {
        Self {
            dist_sq,
            pbc_shift,
            sym_idx,
        }
    }

    /// Return the selected squared Cartesian distance.
    #[must_use]
    pub const fn dist_sq(&self) -> f64 {
        self.dist_sq
    }

    /// Return the selected periodic-cell shift in x, y, z order.
    #[must_use]
    pub const fn pbc_shift(&self) -> &[i32; 3] {
        &self.pbc_shift
    }

    /// Return Gemmi's symmetry-image index; zero denotes the asymmetric unit.
    #[must_use]
    pub const fn sym_idx(&self) -> i32 {
        self.sym_idx
    }
}

// Gemmi❗✔️: struct NearestImage {
// Gemmi❗✔️:   double dist_sq;
// Gemmi❗✔️:   int pbc_shift[3] = { 0, 0, 0 };
// Gemmi❗✔️:   int sym_idx = 0;
// Gemmi❗✔️: };
// Behavior review: this is the typed result shape; only the later nearest-image
// dispatcher defines its initialized observable state. The Rust constructor
// requires all fields rather than representing Gemmi's uninitialized distance.
// Complexity review: the result is three scalar/array fields with no allocation.

fn fractionalize_position(crystal: &BioCrystalInfo, position: [f64; 3]) -> [f64; 3] {
    // Gemmi❗✔️: Fractional fractionalize(const Position& o) const {
    // Gemmi❗✔️:   return Fractional(frac.apply(o));
    // Gemmi❗✔️: }
    // Behavior review: BioTransform::apply retains the full affine fractional
    // translation; unlike a displacement transform, this maps a position.
    // Complexity review: one fixed 3x3 matrix-vector product and translation.
    crystal.fractional().apply(position)
}

fn orthogonalize_difference(crystal: &BioCrystalInfo, fractional_difference: [f64; 3]) -> [f64; 3] {
    // Gemmi❗✔️: Position orthogonalize_difference(const Fractional& delta) const {
    // Gemmi❗✔️:   return Position(orth.mat.multiply(delta));
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Vec3 multiply(const Vec3& p) const {
    // Gemmi❗✔️:   return {a[0][0] * p.x + a[0][1] * p.y + a[0][2] * p.z,
    // Gemmi❗✔️:           a[1][0] * p.x + a[1][1] * p.y + a[1][2] * p.z,
    // Gemmi❗✔️:           a[2][0] * p.x + a[2][1] * p.y + a[2][2] * p.z};
    // Gemmi❗✔️: }
    // Behavior review: multiply the displacement by the orthogonal matrix only;
    // do not apply the transform translation. The canonical BIO matrix-vector
    // helper retains Gemmi's per-row multiplication/addition order.
    // Complexity review: three fixed-length dot products and no allocation.
    multiply_matrix_vector(*crystal.orthogonal().matrix(), fractional_difference)
}

fn squared_length(vector: [f64; 3]) -> f64 {
    // Gemmi❗✔️: Real length_sq() const { return x * x + y * y + z * z; }
    // Behavior review: explicit left-associated x, y, z products and additions.
    // Complexity review: three multiplications and two additions.
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
}

fn iround(value: f64) -> i32 {
    // Gemmi❗✔️: inline int iround(double d) { return static_cast<int>(std::round(d)); }
    // Behavior review: ties round away from zero for finite values whose rounded
    // result fits i32. C++ floating-to-int conversion outside that domain is
    // undefined; Rust's saturating cast there is not claimed as source parity.
    // Complexity review: one rounding operation and one scalar conversion.
    value.round() as i32
}

fn search_pbc_images(
    crystal: &BioCrystalInfo,
    mut difference: [f64; 3],
    image: &mut BioNearestImage,
) -> bool {
    // Gemmi❗✔️: bool search_pbc_images(Fractional&& diff, NearestImage& image) const {
    // Gemmi❗✔️:   int neg_shift[3] = {0, 0, 0};
    // Gemmi❗✔️:   if (is_crystal()) {
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       neg_shift[j] = iround(diff.at(j));
    // Gemmi❗✔️:     diff.x -= neg_shift[0];
    // Gemmi❗✔️:     diff.y -= neg_shift[1];
    // Gemmi❗✔️:     diff.z -= neg_shift[2];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   Position orth_diff = orthogonalize_difference(diff);
    // Gemmi❗✔️:   double dsq = orth_diff.length_sq();
    // Gemmi❗✔️:   if (dsq < image.dist_sq) {
    // Gemmi❗✔️:     image.dist_sq = dsq;
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       image.pbc_shift[j] = -neg_shift[j];
    // Gemmi❗✔️:     return true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return false;
    // Gemmi❗✔️: }
    // Behavior review: this preserves the crystal-only reduction, then scores
    // the residual through orthogonalize_difference (without translation).
    // The image and negated periodic shift change only on a strict improvement.
    // Inputs whose rounded shift cannot be represented as C++ int remain
    // outside the source-defined behavior; iround documents that boundary.
    // Complexity review: fixed three-axis loops, fixed 3x3 arithmetic and no
    // allocation, matching the source's constant-size work and storage.
    let mut negative_shift = [0_i32; 3];
    if crystal.is_crystal() {
        for axis in 0..3 {
            negative_shift[axis] = iround(difference[axis]);
            difference[axis] -= f64::from(negative_shift[axis]);
        }
    }

    let orthogonal_difference = orthogonalize_difference(crystal, difference);
    let distance_sq = squared_length(orthogonal_difference);
    if distance_sq < image.dist_sq {
        image.dist_sq = distance_sq;
        for axis in 0..3 {
            image.pbc_shift[axis] = -negative_shift[axis];
        }
        return true;
    }
    false
}

fn cartesian_distance_sq(reference: [f64; 3], position: [f64; 3]) -> f64 {
    // Gemmi❗✔️: Real dist_sq(const Vec3_& o) const { return (*this - o).length_sq(); }
    // Gemmi❗✔️: Vec3_ operator-(const Vec3_& o) const { return {x-o.x, y-o.y, z-o.z}; }
    // Gemmi❗✔️: Real length_sq() const { return x * x + y * y + z * z; }
    // Behavior review: subtract reference minus position componentwise, then
    // use the same explicit x/y/z squared-length order as the source helpers.
    // Complexity review: three subtractions, three multiplications and two
    // additions; fixed work with no allocation.
    squared_length([
        reference[0] - position[0],
        reference[1] - position[1],
        reference[2] - position[2],
    ])
}

/// Return Gemmi's nearest periodic/symmetry image and its squared distance.
#[must_use]
pub fn find_nearest_image(
    crystal: &BioCrystalInfo,
    reference: [f64; 3],
    position: [f64; 3],
    asu: BioAsu,
) -> BioNearestImage {
    // Gemmi❗✔️: NearestImage find_nearest_image(const Position& ref, const Position& pos, Asu asu) const {
    // Gemmi❗✔️:   NearestImage image;
    // Gemmi❗✔️:   if (asu == Asu::Different)
    // Gemmi❗✔️:     image.dist_sq = INFINITY;
    // Gemmi❗✔️:   else
    // Gemmi❗✔️:     image.dist_sq = ref.dist_sq(pos);
    // Gemmi❗✔️:   if (asu == Asu::Same)
    // Gemmi❗✔️:     return image;
    // Gemmi❗✔️:   Fractional fpos = fractionalize(pos);
    // Gemmi❗✔️:   Fractional fref = fractionalize(ref);
    // Gemmi❗✔️:   search_pbc_images(fpos - fref, image);
    // Gemmi❗✔️:   if (asu == Asu::Different &&
    // Gemmi❗✔️:       image.pbc_shift[0] == 0 && image.pbc_shift[1] == 0 && image.pbc_shift[2] == 0)
    // Gemmi❗✔️:     image.dist_sq = INFINITY;
    // Gemmi❗✔️:   for (int n = 0; n != static_cast<int>(images.size()); ++n)
    // Gemmi❗✔️:     if (search_pbc_images(images[n].apply(fpos) - fref, image))
    // Gemmi❗✔️:       image.sym_idx = n + 1;
    // Gemmi❗✔️:   return image;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Fractional fractionalize(const Position& o) const {
    // Gemmi❗✔️:   return Fractional(frac.apply(o));
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }
    // Gemmi❗✔️: bool search_pbc_images(Fractional&& diff, NearestImage& image) const {
    // Gemmi❗✔️:   int neg_shift[3] = {0, 0, 0};
    // Gemmi❗✔️:   if (is_crystal()) {
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       neg_shift[j] = iround(diff.at(j));
    // Gemmi❗✔️:     diff.x -= neg_shift[0];
    // Gemmi❗✔️:     diff.y -= neg_shift[1];
    // Gemmi❗✔️:     diff.z -= neg_shift[2];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   Position orth_diff = orthogonalize_difference(diff);
    // Gemmi❗✔️:   double dsq = orth_diff.length_sq();
    // Gemmi❗✔️:   if (dsq < image.dist_sq) {
    // Gemmi❗✔️:     image.dist_sq = dsq;
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       image.pbc_shift[j] = -neg_shift[j];
    // Gemmi❗✔️:     return true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return false;
    // Gemmi❗✔️: }
    // Behavior review: retain source initialization, the Same early return,
    // position-fractionalization order, Difference's zero-shift infinity reset,
    // ordered full-image scan and strict candidate replacement. Only a winning
    // symmetry image assigns its one-based index; the source's `images.size()`
    // to `int` loop bound is modeled for source-representable image counts.
    // Complexity review: constant setup plus a single ordered O(image-count)
    // scan; each candidate does fixed-size transform and distance arithmetic,
    // with no image-list clone or allocation.
    let mut image = BioNearestImage::from_parts(
        if asu == BioAsu::Different {
            f64::INFINITY
        } else {
            cartesian_distance_sq(reference, position)
        },
        [0; 3],
        0,
    );
    if asu == BioAsu::Same {
        return image;
    }

    let fractional_position = fractionalize_position(crystal, position);
    let fractional_reference = fractionalize_position(crystal, reference);
    let direct_difference = [
        fractional_position[0] - fractional_reference[0],
        fractional_position[1] - fractional_reference[1],
        fractional_position[2] - fractional_reference[2],
    ];
    let _ = search_pbc_images(crystal, direct_difference, &mut image);
    if asu == BioAsu::Different && image.pbc_shift == [0; 3] {
        image.dist_sq = f64::INFINITY;
    }

    for (index, symmetry_image) in crystal.symmetry_images().iter().enumerate() {
        let transformed_position = symmetry_image.apply(fractional_position);
        let transformed_difference = [
            transformed_position[0] - fractional_reference[0],
            transformed_position[1] - fractional_reference[1],
            transformed_position[2] - fractional_reference[2],
        ];
        if search_pbc_images(crystal, transformed_difference, &mut image) {
            image.sym_idx = index as i32 + 1;
        }
    }
    image
}

#[cfg(test)]
mod tests {
    use super::{
        BioNearestImage, cartesian_distance_sq, find_nearest_image, fractionalize_position, iround,
        orthogonalize_difference, search_pbc_images, squared_length,
    };
    use crate::{BioAsu, BioCrystalCell, BioCrystalInfo, BioTransform};

    fn explicit_crystal(fractional: BioTransform, orthogonal: BioTransform) -> BioCrystalInfo {
        BioCrystalInfo::new(
            BioCrystalCell {
                a: 10.0,
                b: 12.0,
                c: 14.0,
                alpha: 80.0,
                beta: 90.0,
                gamma: 60.0,
            },
            Some("P 1".to_owned()),
            None,
            orthogonal,
            fractional,
            true,
            0,
            Vec::new(),
        )
    }

    fn nearest_image_crystal(symmetry_images: Vec<BioTransform>) -> BioCrystalInfo {
        BioCrystalInfo::new(
            BioCrystalCell {
                a: 8.0,
                b: 8.0,
                c: 8.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
            Some("P 1".to_owned()),
            None,
            BioTransform::new(
                [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[0.125, 0.0, 0.0], [0.0, 0.125, 0.0], [0.0, 0.0, 0.125]],
                [0.0; 3],
            ),
            true,
            0,
            symmetry_images,
        )
    }

    #[test]
    fn lattice_iround_uses_source_half_away_from_zero_rounding() {
        // Gemmi's defined finite, representable domain includes these ties.
        // NaN and rounded values outside i32 are excluded because the pinned
        // C++ floating-to-int conversion is undefined there.
        assert_eq!(iround(0.5), 1);
        assert_eq!(iround(-0.5), -1);
        assert_eq!(iround(1.5), 2);
        assert_eq!(iround(-1.5), -2);
    }

    #[test]
    fn lattice_fractionalization_applies_affine_translation_to_positions() {
        let fractional = BioTransform::new(
            [[2.0, 1.0, 0.0], [0.0, 3.0, -1.0], [1.0, 0.0, 4.0]],
            [7.0, -11.0, 13.0],
        );
        let crystal = explicit_crystal(fractional, BioTransform::identity());

        let actual = fractionalize_position(&crystal, [0.5, -2.0, 1.25]);
        assert_eq!(actual, [6.0, -18.25, 18.5]);
    }

    #[test]
    fn lattice_translated_affine_and_skew_cell_preserve_source_distance_arithmetic() {
        let fractional = BioTransform::new(
            [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
            [7.0, -11.0, 13.0],
        );
        let orthogonal = BioTransform::new(
            [[2.0, 1.0, 0.0], [0.0, 3.0, 0.5], [0.0, 0.0, 4.0]],
            [101.0, -202.0, 303.0],
        );
        let crystal = explicit_crystal(fractional, orthogonal);
        assert!(crystal.is_crystal());
        let reference = [0.0, 0.0, 0.0];
        let position = [1.0, -1.0, 3.0];

        let fractional_reference = fractionalize_position(&crystal, reference);
        let fractional_position = fractionalize_position(&crystal, position);
        assert_eq!(fractional_reference, [7.0, -11.0, 13.0]);
        assert_eq!(fractional_position, [7.5, -11.5, 14.5]);
        assert_eq!(
            cartesian_distance_sq(reference, position).to_bits(),
            11.0_f64.to_bits()
        );

        let difference = [
            fractional_position[0] - fractional_reference[0],
            fractional_position[1] - fractional_reference[1],
            fractional_position[2] - fractional_reference[2],
        ];
        let mut image = BioNearestImage::from_parts(f64::INFINITY, [4, 5, 6], 13);
        assert!(search_pbc_images(&crystal, difference, &mut image));
        assert_eq!(image.pbc_shift(), &[-1, 1, -2]);
        assert_eq!(image.dist_sq().to_bits(), 5.8125_f64.to_bits());
        assert_eq!(image.sym_idx(), 13);
    }

    #[test]
    fn lattice_skew_displacement_omits_translation_and_squares_in_source_order() {
        let orthogonal = BioTransform::new(
            [[2.0, 0.5, 0.0], [0.0, 3.0, 0.25], [0.0, 0.0, 4.0]],
            [101.0, -202.0, 303.0],
        );
        let crystal = explicit_crystal(BioTransform::identity(), orthogonal);
        let fractional_difference = [0.5, -1.5, 2.25];

        let cartesian_difference = orthogonalize_difference(&crystal, fractional_difference);
        assert_eq!(cartesian_difference, [0.25, -3.9375, 9.0]);
        let distance_sq = squared_length(cartesian_difference);
        assert_eq!(distance_sq.to_bits(), 96.56640625_f64.to_bits());

        let result = BioNearestImage::from_parts(distance_sq, [1, -2, 3], 4);
        assert_eq!(result.dist_sq().to_bits(), 96.56640625_f64.to_bits());
        assert_eq!(result.pbc_shift(), &[1, -2, 3]);
        assert_eq!(result.sym_idx(), 4);
    }

    #[test]
    fn lattice_search_pbc_images_reduces_each_crystal_axis_and_negates_shift() {
        let fractional = BioTransform::new(
            [
                [0.1, 0.0, 0.0],
                [0.0, 1.0 / 12.0, 0.0],
                [0.0, 0.0, 1.0 / 14.0],
            ],
            [0.25, -0.5, 1.0],
        );
        let orthogonal = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [101.0, -202.0, 303.0],
        );
        let crystal = explicit_crystal(fractional, orthogonal);
        assert!(crystal.is_crystal());

        for (difference, expected_shift, expected_distance_sq) in [
            ([1.25_f64, 0.0, 0.0], [-1_i32, 0, 0], 0.25_f64),
            ([0.0, -2.25, 0.0], [0, 2, 0], 0.5625_f64),
            ([0.0, 0.0, 3.25], [0, 0, -3], 1.0_f64),
        ] {
            let mut image = BioNearestImage::from_parts(f64::INFINITY, [7, 8, 9], 11);

            assert!(search_pbc_images(&crystal, difference, &mut image));
            assert_eq!(image.pbc_shift(), &expected_shift);
            assert_eq!(image.dist_sq().to_bits(), expected_distance_sq.to_bits());
            assert_eq!(image.sym_idx(), 11);
        }
    }

    #[test]
    fn lattice_search_pbc_images_uses_source_half_ties_on_all_axes() {
        let fractional = BioTransform::new(
            [
                [0.1, 0.0, 0.0],
                [0.0, 1.0 / 12.0, 0.0],
                [0.0, 0.0, 1.0 / 14.0],
            ],
            [0.0; 3],
        );
        let orthogonal = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [0.0; 3],
        );
        let crystal = explicit_crystal(fractional, orthogonal);
        let mut image = BioNearestImage::from_parts(f64::INFINITY, [0; 3], 0);

        assert!(search_pbc_images(&crystal, [0.5, -0.5, 1.5], &mut image));
        assert_eq!(image.pbc_shift(), &[-1, 1, -2]);
        assert_eq!(image.dist_sq().to_bits(), 7.25_f64.to_bits());
    }

    #[test]
    fn lattice_search_pbc_images_does_not_reduce_noncrystal_differences() {
        let noncrystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            None,
            None,
            BioTransform::new(
                [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
                [100.0, 200.0, 300.0],
            ),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        assert!(!noncrystal.is_crystal());
        let mut image = BioNearestImage::from_parts(f64::INFINITY, [4, 5, 6], 7);

        assert!(search_pbc_images(&noncrystal, [1.25, 0.0, 0.0], &mut image));
        assert_eq!(image.pbc_shift(), &[0, 0, 0]);
        assert_eq!(image.dist_sq().to_bits(), 6.25_f64.to_bits());
        assert_eq!(image.sym_idx(), 7);
    }

    #[test]
    fn lattice_search_pbc_images_keeps_existing_result_on_equal_distance() {
        let noncrystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            None,
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        let mut image = BioNearestImage::from_parts(1.0, [9, -8, 7], 42);

        assert!(!search_pbc_images(&noncrystal, [1.0, 0.0, 0.0], &mut image));
        assert_eq!(image.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(image.pbc_shift(), &[9, -8, 7]);
        assert_eq!(image.sym_idx(), 42);
    }

    #[test]
    fn lattice_search_pbc_images_scores_skew_residual_without_translation() {
        let fractional = BioTransform::new(
            [
                [0.1, 0.0, 0.0],
                [0.0, 1.0 / 12.0, 0.0],
                [0.0, 0.0, 1.0 / 14.0],
            ],
            [0.0; 3],
        );
        let orthogonal = BioTransform::new(
            [[2.0, 0.5, 0.0], [0.0, 3.0, 0.25], [0.0, 0.0, 4.0]],
            [101.0, -202.0, 303.0],
        );
        let crystal = explicit_crystal(fractional, orthogonal);
        let mut image = BioNearestImage::from_parts(f64::INFINITY, [0; 3], 99);

        assert!(search_pbc_images(&crystal, [0.5, -1.5, 2.25], &mut image));
        assert_eq!(image.pbc_shift(), &[-1, 2, -2]);
        assert_eq!(image.dist_sq().to_bits(), 4.00390625_f64.to_bits());
        assert_eq!(image.sym_idx(), 99);
    }

    #[test]
    fn lattice_find_nearest_image_preserves_all_asu_branches() {
        let crystal = nearest_image_crystal(Vec::new());
        let reference = [0.0, 0.0, 0.0];
        let adjacent_cell_position = [7.0, 0.0, 0.0];

        let same = find_nearest_image(&crystal, reference, adjacent_cell_position, BioAsu::Same);
        assert_eq!(same.dist_sq().to_bits(), 49.0_f64.to_bits());
        assert_eq!(same.pbc_shift(), &[0, 0, 0]);
        assert_eq!(same.sym_idx(), 0);

        let any = find_nearest_image(&crystal, reference, adjacent_cell_position, BioAsu::Any);
        assert_eq!(any.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(any.pbc_shift(), &[-1, 0, 0]);
        assert_eq!(any.sym_idx(), 0);

        let different = find_nearest_image(
            &crystal,
            reference,
            adjacent_cell_position,
            BioAsu::Different,
        );
        assert_eq!(different.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(different.pbc_shift(), &[-1, 0, 0]);
        assert_eq!(different.sym_idx(), 0);

        let same_cell_different =
            find_nearest_image(&crystal, reference, reference, BioAsu::Different);
        assert_eq!(same_cell_different.dist_sq(), f64::INFINITY);
        assert_eq!(same_cell_different.pbc_shift(), &[0, 0, 0]);
        assert_eq!(same_cell_different.sym_idx(), 0);
    }

    #[test]
    fn lattice_find_nearest_image_uses_the_one_based_index_of_a_later_winner() {
        let inversion = BioTransform::new(
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
            [0.0; 3],
        );
        let inversion_about_quarter = BioTransform::new(
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
            [0.5, 0.0, 0.0],
        );
        let crystal = nearest_image_crystal(vec![inversion, inversion_about_quarter]);

        let image = find_nearest_image(&crystal, [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], BioAsu::Any);
        assert_eq!(image.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(image.pbc_shift(), &[0, 0, 0]);
        assert_eq!(image.sym_idx(), 2);
    }

    #[test]
    fn lattice_find_nearest_image_keeps_the_first_symmetry_candidate_on_a_tie() {
        let inversion_about_quarter = BioTransform::new(
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
            [0.5, 0.0, 0.0],
        );
        let half_cell_translation = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.5, 0.0, 0.0],
        );
        let crystal = nearest_image_crystal(vec![inversion_about_quarter, half_cell_translation]);

        let image = find_nearest_image(&crystal, [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], BioAsu::Any);
        assert_eq!(image.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(image.pbc_shift(), &[0, 0, 0]);
        assert_eq!(image.sym_idx(), 1);
    }
}
