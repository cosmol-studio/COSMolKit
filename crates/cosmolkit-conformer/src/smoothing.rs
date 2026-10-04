//! Dense distance-bound triangle smoothing.

// Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
// @@ All Rights Reserved @@
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

use crate::bounds::BoundsMatrix;

// BEGIN RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds (TriangleSmooth.cpp:14-16)
// RDKit❗✔️: bool triangleSmoothBounds(BoundsMatPtr boundsMat, double tol) {
// RDKit❗✔️:   return triangleSmoothBounds(boundsMat.get(), tol);
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds
pub(super) fn triangle_smooth_bounds_shared(bounds_mat: &mut BoundsMatrix, tol: f64) -> bool {
    triangle_smooth_bounds_ptr(bounds_mat, tol)
}

// BEGIN RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds (TriangleSmooth.cpp:17-70)
// RDKit❗✔️: bool triangleSmoothBounds(BoundsMatrix *boundsMat, double tol) {
// RDKit❗✔️:   auto npt = boundsMat->numRows();
// RDKit❗✔️:   for (auto k = 0u; k < npt; k++) {
// RDKit❗✔️:     for (auto i = 0u; i < npt - 1; i++) {
// RDKit❗✔️:       if (i == k) {
// RDKit❗✔️:         continue;
// RDKit❗✔️:       }
// RDKit❗✔️:       auto ii = i;
// RDKit❗✔️:       auto ik = k;
// RDKit❗✔️:       if (ii > ik) {
// RDKit❗✔️:         std::swap(ii, ik);
// RDKit❗✔️:       }
// RDKit❗✔️:
// RDKit❗✔️:       const auto Uik = boundsMat->getValUnchecked(ii, ik);  // upper bound
// RDKit❗✔️:       const auto Lik = boundsMat->getValUnchecked(ik, ii);  // lower bound
// RDKit❗✔️:       for (auto j = i + 1; j < npt; j++) {
// RDKit❗✔️:         if (j == k) {
// RDKit❗✔️:           continue;
// RDKit❗✔️:         }
// RDKit❗✔️:         auto jj = j;
// RDKit❗✔️:         auto jk = k;
// RDKit❗✔️:         if (jj > jk) {
// RDKit❗✔️:           std::swap(jj, jk);
// RDKit❗✔️:         }
// RDKit❗✔️:         const auto Ukj = boundsMat->getValUnchecked(jj, jk);  // upper bound
// RDKit❗✔️:         const auto sumUikUkj = Uik + Ukj;
// RDKit❗✔️:         if (boundsMat->getValUnchecked(i, j) > sumUikUkj) {
// RDKit❗✔️:           // adjust the upper bound
// RDKit❗✔️:           boundsMat->setValUnchecked(i, j, sumUikUkj);
// RDKit❗✔️:         }
// RDKit❗✔️:
// RDKit❗✔️:         const auto diffLikUjk = Lik - Ukj;
// RDKit❗✔️:         const auto diffLjkUik = boundsMat->getValUnchecked(jk, jj) - Uik;
// RDKit❗✔️:         if (boundsMat->getValUnchecked(j, i) < diffLikUjk) {
// RDKit❗✔️:           // adjust the lower bound
// RDKit❗✔️:           boundsMat->setValUnchecked(j, i, diffLikUjk);
// RDKit❗✔️:         } else if (boundsMat->getValUnchecked(j, i) < diffLjkUik) {
// RDKit❗✔️:           // adjust the lower bound
// RDKit❗✔️:           boundsMat->setValUnchecked(j, i, diffLjkUik);
// RDKit❗✔️:         }
// RDKit❗✔️:         const auto lBound = boundsMat->getValUnchecked(j, i);
// RDKit❗✔️:         const auto uBound = boundsMat->getValUnchecked(i, j);
// RDKit❗✔️:         if (tol > 0. && (lBound - uBound) / lBound > 0. &&
// RDKit❗✔️:             (lBound - uBound) / lBound < tol) {
// RDKit❗✔️:           // adjust the upper bound
// RDKit❗✔️:           boundsMat->setValUnchecked(i, j, lBound);
// RDKit❗✔️:         } else if (lBound - uBound > 0.) {
// RDKit❗✔️:           return false;
// RDKit❗✔️:         }
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️:   return true;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds
//
// The loop index proof and BoundsMatrix's checked square allocation make these
// unchecked accesses safe at each call. They retain RDKit's dense O(1) entry
// cost instead of adding repeated Result and bounds branches inside O(n^3) work.
pub(super) fn triangle_smooth_bounds_ptr(bounds_mat: &mut BoundsMatrix, tol: f64) -> bool {
    let npt = bounds_mat.dimension();
    for k in 0..npt {
        for i in 0..(npt - 1) {
            if i == k {
                continue;
            }
            let (ii, ik) = if i > k { (k, i) } else { (i, k) };

            // SAFETY: k and i are in 0..npt; ordering them preserves that range.
            let uik = unsafe { bounds_mat.get_val_unchecked(ii, ik) };
            // SAFETY: k and i are in 0..npt; ordering them preserves that range.
            let lik = unsafe { bounds_mat.get_val_unchecked(ik, ii) };
            for j in (i + 1)..npt {
                if j == k {
                    continue;
                }
                let (jj, jk) = if j > k { (k, j) } else { (j, k) };
                // SAFETY: j and k are in 0..npt; ordering them preserves that range.
                let ukj = unsafe { bounds_mat.get_val_unchecked(jj, jk) };
                let sum_uik_ukj = uik + ukj;
                // SAFETY: i and j are both in 0..npt and are distinct.
                if unsafe { bounds_mat.get_val_unchecked(i, j) } > sum_uik_ukj {
                    // SAFETY: i and j are both in 0..npt and are distinct.
                    unsafe { bounds_mat.set_val_unchecked(i, j, sum_uik_ukj) };
                }

                let diff_lik_ujk = lik - ukj;
                // SAFETY: j and k are in 0..npt; ordering them preserves that range.
                let diff_ljk_uik = unsafe { bounds_mat.get_val_unchecked(jk, jj) } - uik;
                // SAFETY: i and j are both in 0..npt and are distinct.
                if unsafe { bounds_mat.get_val_unchecked(j, i) } < diff_lik_ujk {
                    // SAFETY: i and j are both in 0..npt and are distinct.
                    unsafe { bounds_mat.set_val_unchecked(j, i, diff_lik_ujk) };
                // SAFETY: i and j are both in 0..npt and are distinct.
                } else if unsafe { bounds_mat.get_val_unchecked(j, i) } < diff_ljk_uik {
                    // SAFETY: i and j are both in 0..npt and are distinct.
                    unsafe { bounds_mat.set_val_unchecked(j, i, diff_ljk_uik) };
                }
                // SAFETY: i and j are both in 0..npt and are distinct.
                let l_bound = unsafe { bounds_mat.get_val_unchecked(j, i) };
                // SAFETY: i and j are both in 0..npt and are distinct.
                let u_bound = unsafe { bounds_mat.get_val_unchecked(i, j) };
                if tol > 0.0
                    && (l_bound - u_bound) / l_bound > 0.0
                    && (l_bound - u_bound) / l_bound < tol
                {
                    // SAFETY: i and j are both in 0..npt and are distinct.
                    unsafe { bounds_mat.set_val_unchecked(i, j, l_bound) };
                } else if l_bound - u_bound > 0.0 {
                    return false;
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use crate::bounds::BoundsMatrix;

    use super::{triangle_smooth_bounds_ptr, triangle_smooth_bounds_shared};

    fn set_pair(bounds: &mut BoundsMatrix, first: usize, second: usize, upper: f64, lower: f64) {
        bounds.set_upper(first, second, upper).unwrap();
        bounds.set_lower(first, second, lower).unwrap();
    }

    #[test]
    fn cf3d_c02_empty_and_one_point_follow_source_loop_guards() {
        let mut empty = BoundsMatrix::new(0).unwrap();
        assert!(triangle_smooth_bounds_ptr(&mut empty, 0.0));

        let mut one_point = BoundsMatrix::new(1).unwrap();
        assert!(triangle_smooth_bounds_shared(&mut one_point, f64::NAN));
    }

    #[test]
    fn cf3d_c02_already_consistent_bounds_remain_unchanged() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        assert!(triangle_smooth_bounds_shared(&mut bounds, 0.25));
        assert_eq!(bounds.get_upper(0, 1).unwrap(), 0.0);
        assert_eq!(bounds.get_lower(0, 1).unwrap(), 0.0);
        assert_eq!(bounds.get_upper(1, 2).unwrap(), 0.0);
        assert_eq!(bounds.get_lower(1, 2).unwrap(), 0.0);
    }

    #[test]
    fn cf3d_c02_upper_bound_tightens_by_strict_sum() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        set_pair(&mut bounds, 0, 1, 2.0, 0.0);
        set_pair(&mut bounds, 0, 2, 3.0, 0.0);
        set_pair(&mut bounds, 1, 2, 20.0, 0.0);

        assert!(triangle_smooth_bounds_ptr(&mut bounds, 0.0));
        assert_eq!(bounds.get_upper(1, 2).unwrap(), 5.0);
        assert_eq!(bounds.get_lower(1, 2).unwrap(), 0.0);
    }

    #[test]
    fn cf3d_c02_first_lower_candidate_wins() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        set_pair(&mut bounds, 0, 1, 8.0, 8.0);
        set_pair(&mut bounds, 0, 2, 2.0, 0.0);
        set_pair(&mut bounds, 1, 2, 20.0, 0.0);

        assert!(triangle_smooth_bounds_shared(&mut bounds, 0.0));
        assert_eq!(bounds.get_upper(1, 2).unwrap(), 10.0);
        assert_eq!(bounds.get_lower(1, 2).unwrap(), 6.0);
    }

    #[test]
    fn cf3d_c02_second_lower_candidate_runs_only_after_first_fails() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        set_pair(&mut bounds, 0, 1, 5.0, 5.0);
        set_pair(&mut bounds, 0, 2, 12.0, 11.0);
        set_pair(&mut bounds, 1, 2, 20.0, 5.0);

        assert!(triangle_smooth_bounds_ptr(&mut bounds, 0.0));
        assert_eq!(bounds.get_upper(1, 2).unwrap(), 17.0);
        assert_eq!(bounds.get_lower(1, 2).unwrap(), 6.0);
    }

    #[test]
    fn cf3d_c02_inconsistency_returns_false_after_prior_upper_write() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        set_pair(&mut bounds, 0, 1, 1.0, 0.0);
        set_pair(&mut bounds, 0, 2, 1.0, 0.0);
        set_pair(&mut bounds, 1, 2, 5.0, 9.0);

        assert!(!triangle_smooth_bounds_ptr(&mut bounds, 0.0));
        assert_eq!(bounds.get_upper(1, 2).unwrap(), 2.0);
        assert_eq!(bounds.get_lower(1, 2).unwrap(), 9.0);
    }

    fn near_threshold_bounds() -> BoundsMatrix {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        set_pair(&mut bounds, 0, 1, 1.0, 0.0);
        set_pair(&mut bounds, 0, 2, 3.0, 0.0);
        set_pair(&mut bounds, 1, 2, 3.0, 3.25);
        bounds
    }

    #[test]
    fn cf3d_c02_tolerance_uses_strict_below_equal_above_comparisons() {
        let mut below = near_threshold_bounds();
        assert!(!triangle_smooth_bounds_ptr(&mut below, 0.07));
        assert_eq!(below.get_upper(1, 2).unwrap(), 3.0);
        assert_eq!(below.get_lower(1, 2).unwrap(), 3.25);

        let mut equal = near_threshold_bounds();
        assert!(!triangle_smooth_bounds_ptr(&mut equal, 1.0 / 13.0));
        assert_eq!(equal.get_upper(1, 2).unwrap(), 3.0);
        assert_eq!(equal.get_lower(1, 2).unwrap(), 3.25);

        let mut above = near_threshold_bounds();
        assert!(triangle_smooth_bounds_shared(&mut above, 0.08));
        assert_eq!(above.get_upper(1, 2).unwrap(), 3.25);
        assert_eq!(above.get_lower(1, 2).unwrap(), 3.25);
    }

    #[test]
    fn cf3d_c02_nonpositive_and_nan_tolerance_do_not_clamp() {
        for tolerance in [0.0, -1.0, f64::NAN] {
            let mut bounds = near_threshold_bounds();
            assert!(!triangle_smooth_bounds_ptr(&mut bounds, tolerance));
            assert_eq!(bounds.get_upper(1, 2).unwrap(), 3.0);
            assert_eq!(bounds.get_lower(1, 2).unwrap(), 3.25);
        }
    }
}
