// Copyright (C)  2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the RDKit BSD license.

const FUNCTOL: f64 = 1.0e-4;
const MOVETOL: f64 = 1.0e-7;
const EPS: f64 = 3.0e-8;
const MAXITS: u32 = 200;
const TOLX: f64 = 4.0 * EPS;
const MAXSTEP: f64 = 100.0;

#[derive(Debug, PartialEq)]
pub(super) struct OptimizerSnapshot {
    pub(super) positions: Vec<f64>,
    pub(super) energy: f64,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) enum OptimizerError<E> {
    Evaluation(E),
    BadTolerance,
    BadDirection,
}

fn update_inverse_hessian(
    dim: usize,
    inv_hessian: &mut [f64],
    d_grad: &mut [f64],
    xi: &[f64],
    hess_d_grad: &mut [f64],
) {
    // BEGIN RDKIT CPP FUNCTION BFGSOpt::minimize inverse-Hessian update
    // RDKit✔️✔️: const double EPS = 3e-8;  //!< Default gradient tolerance in the minimizer
    // RDKit✔️✔️: // compute hessian*dGrad:
    // RDKit✔️✔️: double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   double *ivh = &(invHessian[i * dim]);
    // RDKit✔️✔️:   double &hdgradi = hessDGrad[i];
    // RDKit✔️✔️:   double *dgj = dGrad.data();
    // RDKit✔️✔️:   hdgradi = 0.0;
    // RDKit✔️✔️:   for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
    // RDKit✔️✔️:     hdgradi += *ivh * *dgj;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   fac += dGrad[i] * xi[i];
    // RDKit✔️✔️:   fae += dGrad[i] * hessDGrad[i];
    // RDKit✔️✔️:   sumDGrad += dGrad[i] * dGrad[i];
    // RDKit✔️✔️:   sumXi += xi[i] * xi[i];
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (fac > sqrt(EPS * sumDGrad * sumXi)) {
    // RDKit✔️✔️:   fac = 1.0 / fac;
    // RDKit✔️✔️:   double fad = 1.0 / fae;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     dGrad[i] = fac * xi[i] - fad * hessDGrad[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     unsigned int itab = i * dim;
    // RDKit✔️✔️:     double pxi = fac * xi[i], hdgi = fad * hessDGrad[i],
    // RDKit✔️✔️:            dgi = fae * dGrad[i];
    // RDKit✔️✔️:     double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
    // RDKit✔️✔️:     for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
    // RDKit✔️✔️:       invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
    // RDKit✔️✔️:       invHessian[j * dim + i] = invHessian[itab + j];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BFGSOpt::minimize inverse-Hessian update

    let mut fac = 0.0;
    let mut fae = 0.0;
    let mut sum_d_grad = 0.0;
    let mut sum_xi = 0.0;
    for i in 0..dim {
        hess_d_grad[i] = 0.0;
        for j in 0..dim {
            hess_d_grad[i] += inv_hessian[i * dim + j] * d_grad[j];
        }
        fac += d_grad[i] * xi[i];
        fae += d_grad[i] * hess_d_grad[i];
        sum_d_grad += d_grad[i] * d_grad[i];
        sum_xi += xi[i] * xi[i];
    }
    if fac > (EPS * sum_d_grad * sum_xi).sqrt() {
        fac = 1.0 / fac;
        let fad = 1.0 / fae;
        for i in 0..dim {
            d_grad[i] = fac * xi[i] - fad * hess_d_grad[i];
        }

        for i in 0..dim {
            let pxi = fac * xi[i];
            let hdgi = fad * hess_d_grad[i];
            let dgi = fae * d_grad[i];
            for j in i..dim {
                let upper = i * dim + j;
                inv_hessian[upper] += pxi * xi[j] - hdgi * hess_d_grad[j] + dgi * d_grad[j];
                let updated = inv_hessian[upper];
                inv_hessian[j * dim + i] = updated;
            }
        }
    }
}

pub(super) fn minimize<Context, Energy, Gradient, E>(
    context: &mut Context,
    pos: &mut [f64],
    grad_tol: f64,
    num_iters: &mut u32,
    func_val: &mut f64,
    mut func: Energy,
    mut grad_func: Gradient,
    snapshot_freq: u32,
    mut snapshot_vect: Option<&mut Vec<OptimizerSnapshot>>,
    _func_tol: f64,
    max_its: u32,
) -> Result<i32, OptimizerError<E>>
where
    Energy: FnMut(&mut Context, &mut [f64]) -> Result<f64, E>,
    Gradient: FnMut(&mut Context, &mut [f64], &mut [f64]) -> Result<f64, E>,
{
    // BEGIN RDKIT CPP FUNCTION BFGSOpt::minimize snapshot overload (BFGSOpt.h:184-327)
    // RDKit✔️✔️: const int MAXITS = 200;   //!< Default maximum number of iterations
    // RDKit✔️✔️: const double EPS = 3e-8;  //!< Default gradient tolerance in the minimizer
    // RDKit✔️✔️: const double TOLX = 4. * EPS;  //!< Default direction vector tolerance in the minimizer
    // RDKit✔️✔️: const double MAXSTEP = 100.0;  //!< Default maximum step size in the minimizer
    // RDKit✔️✔️: template <typename EnergyFunctor, typename GradientFunctor>
    // RDKit✔️✔️: int minimize(unsigned int dim, double *pos, double gradTol,
    // RDKit✔️✔️:              unsigned int &numIters, double &funcVal, EnergyFunctor func,
    // RDKit✔️✔️:              GradientFunctor gradFunc, unsigned int snapshotFreq,
    // RDKit✔️✔️:              RDKit::SnapshotVect *snapshotVect, double funcTol = TOLX,
    // RDKit✔️✔️:              unsigned int maxIts = MAXITS) {
    // RDKit✔️✔️:   RDUNUSED_PARAM(funcTol);
    // RDKit✔️✔️:   PRECONDITION(pos, "bad input array");
    // RDKit✔️✔️:   PRECONDITION(gradTol > 0, "bad tolerance");
    // RDKit✔️✔️:   std::vector<double> grad(dim);
    // RDKit✔️✔️:   std::vector<double> dGrad(dim);
    // RDKit✔️✔️:   std::vector<double> hessDGrad(dim);
    // RDKit✔️✔️:   std::vector<double> xi(dim);
    // RDKit✔️✔️:   std::vector<double> invHessian(dim * dim, 0);
    // RDKit✔️✔️:   std::unique_ptr<double[]> newPos(new double[dim]);
    // RDKit✔️✔️:   snapshotFreq = std::min(snapshotFreq, maxIts);
    // RDKit✔️✔️:   double fp = func(pos);
    // RDKit✔️✔️:   gradFunc(pos, grad.data());
    // RDKit✔️✔️:   double sum = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     unsigned int itab = i * dim;
    // RDKit✔️✔️:     invHessian[itab + i] = 1.0;
    // RDKit✔️✔️:     xi[i] = -grad[i];
    // RDKit✔️✔️:     sum += pos[i] * pos[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double maxStep = MAXSTEP * std::max(sqrt(sum), static_cast<double>(dim));
    // RDKit✔️✔️:   for (unsigned int iter = 1; iter <= maxIts; ++iter) {
    // RDKit✔️✔️:     numIters = iter;
    // RDKit✔️✔️:     int status = -1;
    // RDKit✔️✔️:     linearSearch(dim, pos, fp, grad.data(), xi.data(), newPos.get(), funcVal, func,
    // RDKit✔️✔️:                  maxStep, status);
    // RDKit✔️✔️:     CHECK_INVARIANT(status >= 0, "bad direction in linearSearch");
    // RDKit✔️✔️:     fp = funcVal;
    // RDKit✔️✔️:     double test = 0.0;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       xi[i] = newPos[i] - pos[i];
    // RDKit✔️✔️:       pos[i] = newPos[i];
    // RDKit✔️✔️:       double temp = fabs(xi[i]) / std::max(fabs(pos[i]), 1.0);
    // RDKit✔️✔️:       if (temp > test) {
    // RDKit✔️✔️:         test = temp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       dGrad[i] = grad[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (test < TOLX) {
    // RDKit✔️✔️:       if (snapshotVect && snapshotFreq) {
    // RDKit✔️✔️:         RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
    // RDKit✔️✔️:         snapshotVect->push_back(s);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     double gradScale = gradFunc(pos, grad.data());
    // RDKit✔️✔️:     test = 0.0;
    // RDKit✔️✔️:     double term = std::max(funcVal * gradScale, 1.0);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       double temp = fabs(grad[i]) * std::max(fabs(pos[i]), 1.0);
    // RDKit✔️✔️:       test = std::max(test, temp);
    // RDKit✔️✔️:       dGrad[i] = grad[i] - dGrad[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     test /= term;
    // RDKit✔️✔️:     if (test < gradTol) {
    // RDKit✔️✔️:       if (snapshotVect && snapshotFreq) {
    // RDKit✔️✔️:         RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
    // RDKit✔️✔️:         snapshotVect->push_back(s);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // compute hessian*dGrad:
    // RDKit✔️✔️:     double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       double *ivh = &(invHessian[i * dim]);
    // RDKit✔️✔️:       double &hdgradi = hessDGrad[i];
    // RDKit✔️✔️:       double *dgj = dGrad.data();
    // RDKit✔️✔️:       hdgradi = 0.0;
    // RDKit✔️✔️:       for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
    // RDKit✔️✔️:         hdgradi += *ivh * *dgj;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       fac += dGrad[i] * xi[i];
    // RDKit✔️✔️:       fae += dGrad[i] * hessDGrad[i];
    // RDKit✔️✔️:       sumDGrad += dGrad[i] * dGrad[i];
    // RDKit✔️✔️:       sumXi += xi[i] * xi[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (fac > sqrt(EPS * sumDGrad * sumXi)) {
    // RDKit✔️✔️:       fac = 1.0 / fac;
    // RDKit✔️✔️:       double fad = 1.0 / fae;
    // RDKit✔️✔️:       for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:         dGrad[i] = fac * xi[i] - fad * hessDGrad[i];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:         unsigned int itab = i * dim;
    // RDKit✔️✔️:         double pxi = fac * xi[i], hdgi = fad * hessDGrad[i],
    // RDKit✔️✔️:                dgi = fae * dGrad[i];
    // RDKit✔️✔️:         double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
    // RDKit✔️✔️:         for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
    // RDKit✔️✔️:           invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
    // RDKit✔️✔️:           invHessian[j * dim + i] = invHessian[itab + j];
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       unsigned int itab = i * dim;
    // RDKit✔️✔️:       xi[i] = 0.0;
    // RDKit✔️✔️:       double &pxi = xi[i];
    // RDKit✔️✔️:       double *ivh = &(invHessian[itab]);
    // RDKit✔️✔️:       double *gj = grad.data();
    // RDKit✔️✔️:       for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
    // RDKit✔️✔️:         pxi -= *ivh * *gj;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (snapshotVect && snapshotFreq && !(iter % snapshotFreq)) {
    // RDKit✔️✔️:       RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
    // RDKit✔️✔️:       snapshotVect->push_back(s);
    // RDKit✔️✔️:       newPos.reset(new double[dim]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 1;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BFGSOpt::minimize snapshot overload

    let dim = pos.len();
    if !(grad_tol > 0.0) {
        return Err(OptimizerError::BadTolerance);
    }
    let snapshot_freq = snapshot_freq.min(max_its);
    let mut grad = vec![0.0; dim];
    let mut d_grad = vec![0.0; dim];
    let mut hess_d_grad = vec![0.0; dim];
    let mut xi = vec![0.0; dim];
    let mut inv_hessian = vec![0.0; dim * dim];
    let mut new_pos = vec![0.0; dim];

    let mut fp = func(context, pos).map_err(OptimizerError::Evaluation)?;
    grad_func(context, pos, &mut grad).map_err(OptimizerError::Evaluation)?;

    let mut sum = 0.0;
    for i in 0..dim {
        inv_hessian[i * dim + i] = 1.0;
        xi[i] = -grad[i];
        sum += pos[i] * pos[i];
    }
    let sqrt_sum = sum.sqrt();
    let dim_as_f64 = dim as f64;
    let max_step_base = if sqrt_sum < dim_as_f64 {
        dim_as_f64
    } else {
        sqrt_sum
    };
    let max_step = MAXSTEP * max_step_base;

    let mut iter = 1_u32;
    while iter <= max_its {
        *num_iters = iter;

        // RDKit's ForceField energy wrapper is a const callable holding the
        // shared field context; borrowing this callback keeps that context
        // safe without copying a raw owner pointer.
        let status = linear_search(
            context,
            pos,
            fp,
            &grad,
            &mut xi,
            &mut new_pos,
            func_val,
            &mut func,
            max_step,
        )
        .map_err(OptimizerError::Evaluation)?;
        if !(status >= 0) {
            return Err(OptimizerError::BadDirection);
        }
        fp = *func_val;

        let mut test = 0.0;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let abs_pos = pos[i].abs();
            let denominator = if abs_pos < 1.0 { 1.0 } else { abs_pos };
            let temp = xi[i].abs() / denominator;
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i];
        }
        if test < TOLX {
            if let Some(snapshots) = snapshot_vect.as_deref_mut()
                && snapshot_freq != 0
            {
                snapshots.push(OptimizerSnapshot {
                    positions: new_pos,
                    energy: fp,
                });
            }
            return Ok(0);
        }

        let grad_scale = grad_func(context, pos, &mut grad).map_err(OptimizerError::Evaluation)?;
        test = 0.0;
        let func_term = *func_val * grad_scale;
        let term = if func_term < 1.0 { 1.0 } else { func_term };
        for i in 0..dim {
            let abs_pos = pos[i].abs();
            let scale = if abs_pos < 1.0 { 1.0 } else { abs_pos };
            let temp = grad[i].abs() * scale;
            if test < temp {
                test = temp;
            }
            d_grad[i] = grad[i] - d_grad[i];
        }
        test /= term;
        if test < grad_tol {
            if let Some(snapshots) = snapshot_vect.as_deref_mut()
                && snapshot_freq != 0
            {
                snapshots.push(OptimizerSnapshot {
                    positions: new_pos,
                    energy: fp,
                });
            }
            return Ok(0);
        }

        update_inverse_hessian(dim, &mut inv_hessian, &mut d_grad, &xi, &mut hess_d_grad);

        for i in 0..dim {
            xi[i] = 0.0;
            for j in 0..dim {
                xi[i] -= inv_hessian[i * dim + j] * grad[j];
            }
        }
        if let Some(snapshots) = snapshot_vect.as_deref_mut()
            && snapshot_freq != 0
            && iter % snapshot_freq == 0
        {
            let positions = std::mem::take(&mut new_pos);
            snapshots.push(OptimizerSnapshot {
                positions,
                energy: fp,
            });
            new_pos = vec![0.0; dim];
        }
        iter = iter.wrapping_add(1);
    }
    Ok(1)
}

fn linear_search<Context, F, E>(
    context: &mut Context,
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    new_val: &mut f64,
    mut func: F,
    max_step: f64,
) -> Result<i32, E>
where
    F: FnMut(&mut Context, &mut [f64]) -> Result<f64, E>,
{
    // BEGIN RDKIT CPP FUNCTION BFGSOpt::linearSearch (BFGSOpt.h:52-159)
    // RDKit✔️✔️: void linearSearch(unsigned int dim, double *oldPt, double oldVal, double *grad,
    // RDKit✔️✔️:                   double *dir, double *newPt, double &newVal,
    // RDKit✔️✔️:                   EnergyFunctor func, double maxStep, int &resCode) {
    // RDKit✔️✔️:   PRECONDITION(oldPt, "bad input array");
    // RDKit✔️✔️:   PRECONDITION(grad, "bad input array");
    // RDKit✔️✔️:   PRECONDITION(dir, "bad input array");
    // RDKit✔️✔️:   PRECONDITION(newPt, "bad input array");
    // RDKit✔️✔️:   const unsigned int MAX_ITER_LINEAR_SEARCH = 1000;
    // RDKit✔️✔️:   double sum = 0.0, slope = 0.0, test = 0.0, lambda = 0.0;
    // RDKit✔️✔️:   double lambda2 = 0.0, lambdaMin = 0.0, tmpLambda = 0.0, val2 = 0.0;
    // RDKit✔️✔️:   resCode = -1;
    // RDKit✔️✔️:   sum = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     sum += dir[i] * dir[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sum = sqrt(sum);
    // RDKit✔️✔️:   if (sum > maxStep) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       dir[i] *= maxStep / sum;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   slope = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     slope += dir[i] * grad[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (slope >= 0.0) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   test = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     double temp = fabs(dir[i]) / std::max(fabs(oldPt[i]), 1.0);
    // RDKit✔️✔️:     if (temp > test) {
    // RDKit✔️✔️:       test = temp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   lambdaMin = MOVETOL / test;
    // RDKit✔️✔️:   lambda = 1.0;
    // RDKit✔️✔️:   unsigned int it = 0;
    // RDKit✔️✔️:   while (it < MAX_ITER_LINEAR_SEARCH) {
    // RDKit✔️✔️:     if (lambda < lambdaMin) {
    // RDKit✔️✔️:       resCode = 1;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       newPt[i] = oldPt[i] + lambda * dir[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     newVal = func(newPt);
    // RDKit✔️✔️:     if (newVal - oldVal <= FUNCTOL * lambda * slope) {
    // RDKit✔️✔️:       resCode = 0;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (it == 0) {
    // RDKit✔️✔️:       tmpLambda = -slope / (2.0 * (newVal - oldVal - slope));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       double rhs1 = newVal - oldVal - lambda * slope;
    // RDKit✔️✔️:       double rhs2 = val2 - oldVal - lambda2 * slope;
    // RDKit✔️✔️:       double a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) /
    // RDKit✔️✔️:                  (lambda - lambda2);
    // RDKit✔️✔️:       double b = (-lambda2 * rhs1 / (lambda * lambda) +
    // RDKit✔️✔️:                   lambda * rhs2 / (lambda2 * lambda2)) /
    // RDKit✔️✔️:                  (lambda - lambda2);
    // RDKit✔️✔️:       if (a == 0.0) {
    // RDKit✔️✔️:         tmpLambda = -slope / (2.0 * b);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         double disc = b * b - 3 * a * slope;
    // RDKit✔️✔️:         if (disc < 0.0) {
    // RDKit✔️✔️:           tmpLambda = 0.5 * lambda;
    // RDKit✔️✔️:         } else if (b <= 0.0) {
    // RDKit✔️✔️:           tmpLambda = (-b + sqrt(disc)) / (3.0 * a);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           tmpLambda = -slope / (b + sqrt(disc));
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (tmpLambda > 0.5 * lambda) {
    // RDKit✔️✔️:         tmpLambda = 0.5 * lambda;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     lambda2 = lambda;
    // RDKit✔️✔️:     val2 = newVal;
    // RDKit✔️✔️:     lambda = std::max(tmpLambda, 0.1 * lambda);
    // RDKit✔️✔️:     ++it;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     newPt[i] = oldPt[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BFGSOpt::linearSearch

    // Rust slices are non-null by construction; this private caller-owned
    // helper keeps the source's index order and constant auxiliary storage.
    const MAX_ITER_LINEAR_SEARCH: usize = 1000;
    let dim = old_pt.len();
    let mut res_code = -1;
    let mut sum = 0.0;
    for i in 0..dim {
        sum += dir[i] * dir[i];
    }
    sum = sum.sqrt();

    if sum > max_step {
        for value in dir.iter_mut() {
            *value *= max_step / sum;
        }
    }

    let mut slope = 0.0;
    for i in 0..dim {
        slope += dir[i] * grad[i];
    }
    if slope >= 0.0 {
        return Ok(res_code);
    }

    let mut test = 0.0;
    for i in 0..dim {
        let abs_old = old_pt[i].abs();
        // Match std::max(abs_old, 1.0)'s first-argument NaN behavior.
        let denominator = if abs_old < 1.0 { 1.0 } else { abs_old };
        let temp = dir[i].abs() / denominator;
        if temp > test {
            test = temp;
        }
    }

    let lambda_min = MOVETOL / test;
    let mut lambda = 1.0;
    let mut lambda2 = 0.0;
    let mut val2 = 0.0;
    let mut it = 0;
    while it < MAX_ITER_LINEAR_SEARCH {
        if lambda < lambda_min {
            res_code = 1;
            break;
        }
        for i in 0..dim {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        *new_val = func(context, new_pt)?;

        if *new_val - old_val <= FUNCTOL * lambda * slope {
            res_code = 0;
            return Ok(res_code);
        }

        let tmp_lambda = if it == 0 {
            -slope / (2.0 * (*new_val - old_val - slope))
        } else {
            let rhs1 = *new_val - old_val - lambda * slope;
            let rhs2 = val2 - old_val - lambda2 * slope;
            let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
            let b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            let mut cubic_lambda = if a == 0.0 {
                -slope / (2.0 * b)
            } else {
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    0.5 * lambda
                } else if b <= 0.0 {
                    (-b + disc.sqrt()) / (3.0 * a)
                } else {
                    -slope / (b + disc.sqrt())
                }
            };
            if cubic_lambda > 0.5 * lambda {
                cubic_lambda = 0.5 * lambda;
            }
            cubic_lambda
        };

        lambda2 = lambda;
        val2 = *new_val;
        let scaled_lambda = 0.1 * lambda;
        // C++ std::max(a, b) returns a when a < b is false, including NaN.
        lambda = if tmp_lambda < scaled_lambda {
            scaled_lambda
        } else {
            tmp_lambda
        };
        it += 1;
    }

    for i in 0..dim {
        new_pt[i] = old_pt[i];
    }
    Ok(res_code)
}

#[cfg(test)]
mod tests {
    use std::convert::Infallible;

    use super::{
        EPS, OptimizerError, OptimizerSnapshot, TOLX, linear_search, minimize,
        update_inverse_hessian,
    };

    fn ok_f64(value: f64) -> Result<f64, Infallible> {
        Ok(value)
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum CallbackFailure {
        InitialEnergy,
        InitialGradient,
        FirstTrialEnergy,
        BacktrackingEnergy,
        IterationGradient,
    }

    #[derive(Default)]
    struct CallbackTrace {
        events: Vec<&'static str>,
        energy_calls: usize,
        gradient_calls: usize,
        mutations: usize,
        observed_trial: Option<f64>,
    }

    impl CallbackTrace {
        fn record(&mut self, event: &'static str) {
            self.events.push(event);
            self.mutations += 1;
        }
    }

    #[test]
    fn cf3d_opt_bridge_invariant_zero_negative_and_nan_tolerance_are_typed() {
        for grad_tol in [0.0, -1.0, f64::NAN] {
            let mut callback_calls = 0;
            let mut positions = [2.0];
            let mut iterations = 17;
            let mut final_energy = 91.0;
            let mut snapshots = vec![OptimizerSnapshot {
                positions: vec![-1.0],
                energy: -1.0,
            }];
            let mut energy = |calls: &mut usize, _: &mut [f64]| -> Result<f64, Infallible> {
                *calls += 1;
                Ok(4.0)
            };
            let mut gradient =
                |calls: &mut usize, _: &mut [f64], grad: &mut [f64]| -> Result<f64, Infallible> {
                    *calls += 1;
                    grad.fill(0.0);
                    Ok(1.0)
                };

            let result = minimize(
                &mut callback_calls,
                &mut positions,
                grad_tol,
                &mut iterations,
                &mut final_energy,
                &mut energy,
                &mut gradient,
                1,
                Some(&mut snapshots),
                0.0,
                4,
            );

            assert_eq!(result, Err(OptimizerError::BadTolerance));
            assert_eq!(callback_calls, 0);
            assert_eq!(positions, [2.0]);
            assert_eq!(iterations, 17);
            assert_eq!(final_energy, 91.0);
            assert_eq!(snapshots.len(), 1);
            assert_eq!(snapshots[0].positions, [-1.0]);
            assert_eq!(snapshots[0].energy, -1.0);
        }
    }

    #[test]
    fn cf3d_opt_bridge_invariant_zero_gradient_is_bad_direction() {
        let mut callback_calls = [0, 0];
        let mut positions = [2.0];
        let mut iterations = 0;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |calls: &mut [usize; 2], _: &mut [f64]| -> Result<f64, Infallible> {
            calls[0] += 1;
            Ok(4.0)
        };
        let mut gradient = |calls: &mut [usize; 2], _: &mut [f64], grad: &mut [f64]| {
            calls[1] += 1;
            grad[0] = 0.0;
            Ok::<f64, Infallible>(1.0)
        };

        let result = minimize(
            &mut callback_calls,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(result, Err(OptimizerError::BadDirection));
        assert_eq!(callback_calls, [1, 1]);
        assert_eq!(positions, [2.0]);
        assert_eq!(iterations, 1);
        assert_eq!(final_energy, 91.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_invariant_callback_error_stays_distinct() {
        let mut trace = CallbackTrace::default();
        let mut positions = [2.0];
        let mut iterations = 17;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, _: &mut [f64]| {
            trace.record("energy.initial");
            Err(CallbackFailure::InitialEnergy)
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], _: &mut [f64]| {
            trace.record("gradient.initial");
            Ok(1.0)
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(CallbackFailure::InitialEnergy))
        );
        assert_ne!(result, Err(OptimizerError::BadTolerance));
        assert_ne!(result, Err(OptimizerError::BadDirection));
        assert_eq!(trace.events, ["energy.initial"]);
        assert_eq!(positions, [2.0]);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_invariant_nonconvergence_keeps_status_and_snapshot_schedule() {
        let mut positions = [2.0];
        let mut iterations = 0;
        let mut final_energy = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            ok_f64(-point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = -300.0;
            ok_f64(1.0)
        };
        let mut snapshots = Vec::new();

        let result = minimize(
            &mut (),
            &mut positions,
            0.1,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            1,
        );

        assert_eq!(result, Ok(1));
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 2);
        assert_eq!(iterations, 1);
        assert_eq!(final_energy, -202.0);
        assert_eq!(positions, [202.0]);
        assert_eq!(snapshots.len(), 1);
        assert_eq!(snapshots[0].positions, [202.0]);
        assert_eq!(snapshots[0].energy, -202.0);
    }

    #[test]
    fn cf3d_opt_bridge_failure_initial_energy_preserves_callback_state() {
        let mut trace = CallbackTrace::default();
        let mut positions = [2.0];
        let mut iterations = 17;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, _: &mut [f64]| {
            trace.record("energy.initial");
            Err(CallbackFailure::InitialEnergy)
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], grad: &mut [f64]| {
            trace.record("gradient.initial");
            grad.fill(0.0);
            Ok(1.0)
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(CallbackFailure::InitialEnergy))
        );
        assert_eq!(trace.events, ["energy.initial"]);
        assert_eq!(trace.mutations, 1);
        assert_eq!(positions, [2.0]);
        assert_eq!(iterations, 17);
        assert_eq!(final_energy, 91.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_failure_initial_gradient_preserves_callback_state() {
        let mut trace = CallbackTrace::default();
        let mut positions = [2.0];
        let mut iterations = 17;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, _: &mut [f64]| {
            trace.record("energy.initial");
            Ok(4.0)
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], grad: &mut [f64]| {
            trace.record("gradient.initial");
            grad[0] = 1.0;
            Err(CallbackFailure::InitialGradient)
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(CallbackFailure::InitialGradient))
        );
        assert_eq!(trace.events, ["energy.initial", "gradient.initial"]);
        assert_eq!(trace.mutations, 2);
        assert_eq!(positions, [2.0]);
        assert_eq!(iterations, 17);
        assert_eq!(final_energy, 91.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_failure_first_trial_energy_stops_at_source_call() {
        let mut trace = CallbackTrace::default();
        let mut positions = [2.0];
        let mut iterations = 17;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, point: &mut [f64]| {
            let call = trace.energy_calls;
            trace.energy_calls += 1;
            if call == 0 {
                trace.record("energy.initial");
                Ok(point[0] * point[0])
            } else {
                trace.record("energy.trial.first");
                trace.observed_trial = Some(point[0]);
                Err(CallbackFailure::FirstTrialEnergy)
            }
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], grad: &mut [f64]| {
            trace.gradient_calls += 1;
            trace.record("gradient.initial");
            grad[0] = 1.0;
            Ok(1.0)
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(
                CallbackFailure::FirstTrialEnergy
            ))
        );
        assert_eq!(
            trace.events,
            ["energy.initial", "gradient.initial", "energy.trial.first"]
        );
        assert_eq!(trace.energy_calls, 2);
        assert_eq!(trace.gradient_calls, 1);
        assert_eq!(trace.mutations, 3);
        assert_eq!(trace.observed_trial, Some(1.0));
        assert_eq!(positions, [2.0]);
        assert_eq!(iterations, 1);
        assert_eq!(final_energy, 91.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_failure_backtracking_energy_preserves_trial_prefix() {
        let mut trace = CallbackTrace::default();
        let mut positions = [0.0];
        let mut iterations = 17;
        let mut final_energy = 91.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, point: &mut [f64]| {
            let call = trace.energy_calls;
            trace.energy_calls += 1;
            match call {
                0 => {
                    trace.record("energy.initial");
                    Ok(0.0)
                }
                1 => {
                    trace.record("energy.trial.first");
                    Ok(1.0)
                }
                _ => {
                    trace.record("energy.trial.backtrack");
                    trace.observed_trial = Some(point[0]);
                    Err(CallbackFailure::BacktrackingEnergy)
                }
            }
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], grad: &mut [f64]| {
            trace.gradient_calls += 1;
            trace.record("gradient.initial");
            grad[0] = 1.0;
            Ok(1.0)
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(
                CallbackFailure::BacktrackingEnergy
            ))
        );
        assert_eq!(
            trace.events,
            [
                "energy.initial",
                "gradient.initial",
                "energy.trial.first",
                "energy.trial.backtrack"
            ]
        );
        assert_eq!(trace.energy_calls, 3);
        assert_eq!(trace.gradient_calls, 1);
        assert_eq!(trace.mutations, 4);
        assert_eq!(trace.observed_trial, Some(-0.25));
        assert_eq!(positions, [0.0]);
        assert_eq!(iterations, 1);
        assert_eq!(final_energy, 1.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_opt_bridge_failure_later_gradient_keeps_prior_snapshot_and_position() {
        let mut trace = CallbackTrace::default();
        let mut positions = [0.0];
        let mut iterations = 0;
        let mut final_energy = 0.0;
        let mut snapshots = Vec::new();
        let mut energy = |trace: &mut CallbackTrace, point: &mut [f64]| {
            let call = trace.energy_calls;
            trace.energy_calls += 1;
            if call == 0 {
                trace.record("energy.initial");
                Ok(-point[0])
            } else {
                trace.record("energy.trial");
                Ok(-point[0])
            }
        };
        let mut gradient = |trace: &mut CallbackTrace, _: &mut [f64], grad: &mut [f64]| {
            let call = trace.gradient_calls;
            trace.gradient_calls += 1;
            match call {
                0 => trace.record("gradient.initial"),
                1 => trace.record("gradient.iteration.first"),
                _ => trace.record("gradient.iteration.second"),
            }
            grad[0] = -1.0;
            if call == 2 {
                Err(CallbackFailure::IterationGradient)
            } else {
                Ok(1.0)
            }
        };

        let result = minimize(
            &mut trace,
            &mut positions,
            0.1,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            1,
            Some(&mut snapshots),
            0.0,
            4,
        );

        assert_eq!(
            result,
            Err(OptimizerError::Evaluation(
                CallbackFailure::IterationGradient
            ))
        );
        assert_eq!(
            trace.events,
            [
                "energy.initial",
                "gradient.initial",
                "energy.trial",
                "gradient.iteration.first",
                "energy.trial",
                "gradient.iteration.second"
            ]
        );
        assert_eq!(trace.energy_calls, 3);
        assert_eq!(trace.gradient_calls, 3);
        assert_eq!(trace.mutations, 6);
        assert_eq!(positions, [2.0]);
        assert_eq!(iterations, 2);
        assert_eq!(final_energy, -2.0);
        assert_eq!(snapshots.len(), 1);
        assert_eq!(snapshots[0].positions, [1.0]);
        assert_eq!(snapshots[0].energy, -1.0);
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "actual={actual:?}, expected={expected:?}"
        );
    }

    // Fixed inverse-Hessian values are derived from RDKit BFGSOpt.h:286-307,
    // pinned file SHA-256 5b6df4743cb79d2515f2d4b3a0c11723b02ef7ac3ec5ab09b3823778db4cc801.
    #[test]
    fn cf3d_f11_guard_not_taken_below_threshold_preserves_matrix_and_delta() {
        let mut inv_hessian = [1.0, 0.0, 0.0, 1.0];
        let mut d_grad = [1.0, 0.0];
        let xi = [0.0, 1.0];
        let mut hess_d_grad = [9.0, 9.0];

        update_inverse_hessian(2, &mut inv_hessian, &mut d_grad, &xi, &mut hess_d_grad);

        assert_eq!(inv_hessian, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(d_grad, [1.0, 0.0]);
        assert_eq!(hess_d_grad, [1.0, 0.0]);
    }

    #[test]
    fn cf3d_f11_exact_positive_threshold_skips_update() {
        let mut inv_hessian = [1.0, 0.0, 0.0, 1.0];
        let mut d_grad = [1.0, 0.0];
        let xi = [EPS.sqrt(), (1.0 - EPS).sqrt()];
        let mut hess_d_grad = [0.0, 0.0];

        let sum_xi = xi[0] * xi[0] + xi[1] * xi[1];
        let fac = d_grad[0] * xi[0] + d_grad[1] * xi[1];
        let threshold = (EPS * sum_xi).sqrt();
        assert_eq!(sum_xi, 1.0);
        assert_eq!(fac, threshold);

        update_inverse_hessian(2, &mut inv_hessian, &mut d_grad, &xi, &mut hess_d_grad);

        assert_eq!(inv_hessian, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(d_grad, [1.0, 0.0]);
        assert_eq!(hess_d_grad, [1.0, 0.0]);
    }

    #[test]
    fn cf3d_f11_guard_taken_matches_reference_matrix_and_enforces_symmetry() {
        let mut inv_hessian = [1.0, 0.0, 0.5, 1.0];
        let mut d_grad = [1.0, 1.0];
        let xi = [1.0, 2.0];
        let mut hess_d_grad = [0.0, 0.0];

        update_inverse_hessian(2, &mut inv_hessian, &mut d_grad, &xi, &mut hess_d_grad);

        assert_eq!(hess_d_grad, [1.0, 1.5]);
        assert_close(d_grad[0], -1.0 / 15.0);
        assert_close(d_grad[1], 1.0 / 15.0);
        assert_close(inv_hessian[0], 17.0 / 18.0);
        assert_close(inv_hessian[1], 1.0 / 18.0);
        assert_close(inv_hessian[2], 1.0 / 18.0);
        assert_close(inv_hessian[3], 13.0 / 9.0);
    }

    #[test]
    fn cf3d_f11_source_fae_zero_reciprocal_has_no_extra_guard() {
        let mut inv_hessian = [0.0];
        let mut d_grad = [1.0];
        let xi = [1.0];
        let mut hess_d_grad = [8.0];

        update_inverse_hessian(1, &mut inv_hessian, &mut d_grad, &xi, &mut hess_d_grad);

        assert!(d_grad[0].is_nan());
        assert!(inv_hessian[0].is_nan());
        assert_eq!(hess_d_grad, [0.0]);
    }

    #[test]
    fn cf3d_f10_first_trial_accepts_source_ordered_point_and_value() {
        // RDKit source: BFGSOpt.h:52-159; quadratic energy is a fixed oracle.
        let old_pt = [1.0];
        let grad = [1.0];
        let mut dir = [-1.0];
        let mut new_pt = [9.0];
        let mut new_val = 17.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            calls += 1;
            ok_f64(point[0] * point[0])
        };

        let status = linear_search(
            &mut (),
            &old_pt,
            1.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(calls, 1);
        assert_eq!(new_pt, [0.0]);
        assert_eq!(new_val, 0.0);
        assert_eq!(dir, [-1.0]);
    }

    #[test]
    fn cf3d_f10_quadratic_backtracking_preserves_floor_and_unclamped_first_step() {
        // RDKit source: BFGSOpt.h:112-133 and 149-152.
        let old_pt = [0.0];
        let grad = [1.0];

        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(if calls == 1 { 1.0 } else { -0.1 })
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 2);
        assert_eq!(new_pt, [-0.25]);
        assert_eq!(new_val, -0.1);

        // The first quadratic proposal is not subject to the cubic half-step
        // cap: raw lambda is 0.500025..., and the next trial accepts it.
        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(if calls == 1 { -0.00005 } else { -1.0 })
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 2);
        assert_close(new_pt[0], -0.5000250012500625);
        assert_eq!(new_val, -1.0);

        // A quadratic proposal below 0.1*lambda uses the source floor.
        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(if calls == 1 { 9.0 } else { -1.0 })
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 2);
        assert_eq!(new_pt, [-0.1]);
        assert_eq!(new_val, -1.0);
    }

    #[test]
    fn cf3d_f10_cubic_backtracking_covers_source_coefficient_branches() {
        // RDKit source: BFGSOpt.h:118-146 and 149-152.
        let old_pt = [0.0];
        let grad = [1.0];

        // With the first two rejected values at zero, a=-2, b=3, and
        // disc=3; the positive-b cubic root is 1/(3+sqrt(3)).
        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let values = [0.0, 0.0, -1.0];
        let mut energy = |_: &mut (), _: &mut [f64]| {
            let value = values[calls];
            calls += 1;
            ok_f64(value)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 3);
        assert_close(new_pt[0], -1.0 / (3.0 + 3.0_f64.sqrt()));
        assert_eq!(new_val, -1.0);

        // These fixed values make a exactly zero in the source expression;
        // the a==0 quadratic root is then subject to the 0.1 floor.
        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let values = [9.0, 1.3877787807814457e-17, -1.0];
        let mut energy = |_: &mut (), _: &mut [f64]| {
            let value = values[calls];
            calls += 1;
            ok_f64(value)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 3);
        assert_eq!(new_pt, [-0.05]);
        assert_eq!(new_val, -1.0);

        // A large first rejection drives b<=0; the source root exceeds
        // 0.5*lambda and is capped before the final accepted trial.
        let mut dir = [-1.0];
        let mut new_pt = [8.0];
        let mut new_val = 22.0;
        let mut calls = 0;
        let values = [999.0, -0.000009, -1.0];
        let mut energy = |_: &mut (), _: &mut [f64]| {
            let value = values[calls];
            calls += 1;
            ok_f64(value)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 0);
        assert_eq!(calls, 3);
        assert_eq!(new_pt, [-0.05]);
        assert_eq!(new_val, -1.0);
    }

    #[test]
    fn cf3d_f10_minimum_lambda_restores_coordinates_and_preserves_values() {
        // RDKit source: BFGSOpt.h:83-108 and 148-159.
        let old_pt = [1.0e308];
        let grad = [1.0];
        let mut dir = [-f64::MIN_POSITIVE];
        let mut new_pt = [7.0];
        let mut new_val = 42.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(0.0)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 1);
        assert_eq!(calls, 0);
        assert_eq!(new_pt, old_pt);
        assert_eq!(new_val, 42.0);

        // A rejected trial proposes a scale below the 0.1 floor; the floored
        // next lambda is below lambdaMin, so status 1 restores the old point
        // while retaining the last evaluated energy.
        let old_pt = [0.0];
        let grad = [1.0];
        let mut dir = [-5.0e-7];
        let mut new_pt = [7.0];
        let mut new_val = 42.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(1.0)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 1);
        assert_eq!(calls, 1);
        assert_eq!(new_pt, old_pt);
        assert_eq!(new_val, 1.0);
    }

    #[test]
    fn cf3d_f10_positive_and_zero_slope_leave_outputs_untouched() {
        // RDKit source: BFGSOpt.h:77-81.
        let old_pt = [2.0];
        let grad = [1.0];
        let mut new_pt = [8.0];
        let mut new_val = 17.0;
        let mut forbidden_energy = |_: &mut (), _: &mut [f64]| -> Result<f64, Infallible> {
            panic!("non-descent directions must return before energy evaluation")
        };

        let mut positive_dir = [1.0];
        assert_eq!(
            linear_search(
                &mut (),
                &old_pt,
                4.0,
                &grad,
                &mut positive_dir,
                &mut new_pt,
                &mut new_val,
                &mut forbidden_energy,
                10.0,
            )
            .unwrap(),
            -1
        );
        assert_eq!(new_pt, [8.0]);
        assert_eq!(new_val, 17.0);

        let mut zero_dir = [0.0];
        assert_eq!(
            linear_search(
                &mut (),
                &old_pt,
                4.0,
                &grad,
                &mut zero_dir,
                &mut new_pt,
                &mut new_val,
                &mut forbidden_energy,
                10.0,
            )
            .unwrap(),
            -1
        );
        assert_eq!(new_pt, [8.0]);
        assert_eq!(new_val, 17.0);
    }

    #[test]
    fn cf3d_f10_max_step_scales_the_mutated_direction_before_slope() {
        // RDKit source: BFGSOpt.h:68-81 and 109-115.
        let old_pt = [0.0, 0.0];
        let grad = [-1.0, -1.0];
        let mut dir = [3.0, 4.0];
        let mut new_pt = [9.0, 9.0];
        let mut new_val = 22.0;
        let mut energy = |_: &mut (), _: &mut [f64]| ok_f64(-1.0);

        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            2.5,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(dir, [1.5, 2.0]);
        assert_eq!(new_pt, [1.5, 2.0]);
        assert_eq!(new_val, -1.0);
    }

    #[test]
    fn cf3d_f10_std_max_nan_order_and_iteration_limit_match_source() {
        // RDKit source: BFGSOpt.h:83-108 and 149-159.
        let old_pt = [f64::NAN];
        let grad = [1.0];
        let mut dir = [-1.0];
        let mut new_pt = [7.0];
        let mut new_val = 42.0;
        let mut calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            calls += 1;
            ok_f64(0.0)
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, 1);
        assert_eq!(calls, 0);
        assert!(new_pt[0].is_nan());
        assert_eq!(new_val, 42.0);

        let old_pt = [0.0];
        let grad = [1.0];
        let mut dir = [-1.0];
        let mut new_pt = [7.0];
        let mut new_val = 42.0;
        let mut calls = 0;
        let mut saw_nan_trial = false;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            calls += 1;
            if calls == 1 {
                ok_f64(f64::NAN)
            } else {
                saw_nan_trial |= point[0].is_nan();
                ok_f64(0.0)
            }
        };
        let status = linear_search(
            &mut (),
            &old_pt,
            0.0,
            &grad,
            &mut dir,
            &mut new_pt,
            &mut new_val,
            &mut energy,
            10.0,
        )
        .unwrap();
        assert_eq!(status, -1);
        assert_eq!(calls, 1000);
        assert!(saw_nan_trial);
        assert_eq!(new_pt, old_pt);
        assert_eq!(new_val, 0.0);
    }

    // Fixed driver expectations follow RDKit BFGSOpt.h:184-327, pinned SHA-256
    // 5b6df4743cb79d2515f2d4b3a0c11723b02ef7ac3ec5ab09b3823778db4cc801.
    #[test]
    fn cf3d_f12_zero_iterations_evaluates_initial_state_and_preserves_outputs() {
        let mut pos = [2.0];
        let mut num_iters = 31;
        let mut func_val = -17.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            ok_f64(point[0] * point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = 4.0;
            ok_f64(2.0)
        };
        let mut snapshots = Vec::new();

        let status = minimize(
            &mut (),
            &mut pos,
            1.0e-6,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            9,
            Some(&mut snapshots),
            123.0,
            0,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(energy_calls, 1);
        assert_eq!(gradient_calls, 1);
        assert_eq!(num_iters, 31);
        assert_eq!(func_val, -17.0);
        assert_eq!(pos, [2.0]);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_f12_position_convergence_records_terminal_snapshot_off_frequency() {
        // This displacement exceeds F10 MOVETOL (1e-7) but stays below F12 TOLX (1.2e-7).
        let mut pos = [1.1e-7];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            ok_f64(0.5 * point[0] * point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), point: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = point[0];
            ok_f64(1.0)
        };
        let mut snapshots = Vec::new();

        let status = minimize(
            &mut (),
            &mut pos,
            1.0e-6,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            2,
            Some(&mut snapshots),
            0.0,
            4,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 1);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, 0.0);
        assert_eq!(pos, [0.0]);
        assert_eq!(snapshots.len(), 1);
        assert_eq!(snapshots[0].positions, [0.0]);
        assert_eq!(snapshots[0].energy, 0.0);
    }

    #[test]
    fn cf3d_f12_gradient_convergence_records_terminal_snapshot() {
        let mut pos = [1.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(point[0]);
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { 1.0 } else { 0.0 };
            ok_f64(1.0)
        };
        let mut snapshots = Vec::new();

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            3,
            Some(&mut snapshots),
            0.0,
            3,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, 0.0);
        assert_eq!(pos, [0.0]);
        assert_eq!(snapshots.len(), 1);
        assert_eq!(snapshots[0].positions, [0.0]);
        assert_eq!(snapshots[0].energy, 0.0);
    }

    #[test]
    fn cf3d_f12_one_iteration_cap_returns_nonconverged_status_without_snapshot_target() {
        let mut pos = [2.0];
        let mut num_iters = 0;
        let mut func_val = -1.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(-point[0]);
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            grad[0] = -300.0;
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            1,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, -202.0);
        assert_eq!(pos, [202.0]);
    }

    #[test]
    fn cf3d_f12_dimension_floor_limits_step_to_one_hundred_times_dimension() {
        let mut pos = [0.5];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), _: &mut [f64]| {
            energy_calls += 1;
            ok_f64(if energy_calls == 1 { 0.0 } else { -1000.0 })
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { -200.0 } else { 0.0 };
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, -1000.0);
        assert_eq!(pos, [100.5]);
    }

    #[test]
    fn cf3d_f12_large_coordinate_position_scaling_converges_before_gradient() {
        let mut pos = [1000.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            ok_f64(-1.1e-4 * point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = -1.1e-4;
            ok_f64(1.0)
        };
        let mut snapshots = Vec::new();

        let status = minimize(
            &mut (),
            &mut pos,
            1.0e-6,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            2,
            Some(&mut snapshots),
            0.0,
            4,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 1);
        assert_eq!(num_iters, 1);
        assert_close(pos[0], 1000.00011);
        assert_eq!(snapshots.len(), 1);
        assert_close(snapshots[0].positions[0], pos[0]);
    }

    #[test]
    fn cf3d_f12_gradient_term_uses_energy_times_gradient_scale_above_one() {
        let mut pos = [0.5];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(4.0 + point[0]);
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { 1.0 } else { 3.5 };
            ok_f64(if gradient_calls == 1 { 1.0 } else { 2.0 })
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.75,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, 3.5);
        assert_eq!(pos, [-0.5]);
    }

    #[test]
    fn cf3d_f12_gradient_test_scales_large_coordinates_before_threshold() {
        let mut pos = [10.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(100.0 + point[0]);
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { 1.0 } else { 0.01 };
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.0005,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, 109.0);
        assert_eq!(pos, [9.0]);
    }

    #[test]
    fn cf3d_f12_zero_snapshot_frequency_suppresses_terminal_snapshot() {
        let mut pos = [1.1e-7];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(0.5 * point[0] * point[0]);
        let mut gradient = |_: &mut (), point: &mut [f64], grad: &mut [f64]| {
            grad[0] = point[0];
            ok_f64(1.0)
        };
        let mut snapshots = Vec::new();

        let status = minimize(
            &mut (),
            &mut pos,
            1.0e-6,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            Some(&mut snapshots),
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(num_iters, 1);
        assert_eq!(pos, [0.0]);
        assert_eq!(func_val, 0.0);
        assert!(snapshots.is_empty());
    }

    #[test]
    fn cf3d_f12_snapshot_frequency_clamps_to_cap_and_appends_after_update() {
        let mut pos = [2.0];
        let mut num_iters = 0;
        let mut func_val = -1.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(point[0]);
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            grad[0] = 1.0;
            ok_f64(1.0)
        };
        let mut snapshots = vec![OptimizerSnapshot {
            positions: vec![-1.0],
            energy: -1.0,
        }];

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            9,
            Some(&mut snapshots),
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(num_iters, 1);
        assert_eq!(snapshots.len(), 2);
        assert_eq!(snapshots[0].positions, [-1.0]);
        assert_eq!(snapshots[0].energy, -1.0);
        assert_eq!(snapshots[1].positions, [1.0]);
        assert_eq!(snapshots[1].energy, 1.0);
    }

    #[test]
    fn cf3d_f12_position_threshold_is_strict() {
        let mut pos = [TOLX];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(0.5 * point[0] * point[0]);
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { TOLX } else { 1.0 };
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.5,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(pos, [0.0]);
    }

    #[test]
    fn cf3d_f12_gradient_threshold_is_strict() {
        let mut pos = [1.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy = |_: &mut (), point: &mut [f64]| ok_f64(point[0]);
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            grad[0] = if gradient_calls == 1 { 1.0 } else { 0.5 };
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.5,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(pos, [0.0]);
    }

    #[test]
    fn cf3d_f12_max_step_nan_order_matches_source() {
        let mut pos = [f64::NAN, 0.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy_calls = 0;
        let mut observed_second_coordinate = f64::NAN;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            if energy_calls == 1 {
                ok_f64(0.0)
            } else {
                observed_second_coordinate = point[1];
                ok_f64(-1000.0)
            }
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            if gradient_calls == 1 {
                grad.copy_from_slice(&[-300.0, -100.0]);
            } else {
                grad.fill(0.0);
            }
            ok_f64(1.0)
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 2);
        assert_eq!(observed_second_coordinate, 100.0);
        assert!(pos[0].is_nan());
        assert_eq!(pos[1], 100.0);
    }

    #[test]
    fn cf3d_f12_gradient_term_nan_order_does_not_converge() {
        let mut pos = [1.0];
        let mut num_iters = 0;
        let mut func_val = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| {
            energy_calls += 1;
            ok_f64(point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (), _: &mut [f64], grad: &mut [f64]| {
            gradient_calls += 1;
            if gradient_calls == 1 {
                grad[0] = 1.0;
                ok_f64(1.0)
            } else {
                grad[0] = 0.0;
                ok_f64(f64::NAN)
            }
        };

        let status = minimize(
            &mut (),
            &mut pos,
            0.1,
            &mut num_iters,
            &mut func_val,
            &mut energy,
            &mut gradient,
            0,
            None,
            0.0,
            1,
        )
        .unwrap();

        assert_eq!(status, 1);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 2);
        assert_eq!(num_iters, 1);
        assert_eq!(func_val, 0.0);
    }
}
