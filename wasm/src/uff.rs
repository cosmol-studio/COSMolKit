//! Thin UFF evaluation and detached optimization projections.
use crate::Molecule;
use cosmolkit as ck;
use std::cell::RefCell;
pub struct UffOptimizationResult {
    inner: ck::UffOptimizationResult,
}
impl UffOptimizationResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule().clone(),
        Molecule {
            inner: RefCell::new(self.inner.molecule().clone()),
        }
    }
    pub fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    pub fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
}
pub struct UffConformerOptimizationResult {
    inner: ck::UffConformerOptimizationResult,
}
impl UffConformerOptimizationResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule().clone(),
        Molecule {
            inner: RefCell::new(self.inner.molecule().clone()),
        }
    }
    pub fn conformer_results(&self) -> &[ck::UffConformerResult] {
        self.inner.conformer_results()
    }
}
impl Molecule {
    pub fn uff_energy_gradient(&self) -> Result<ck::UffEnergyGradient, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.uff_energy_gradient()
        self.inner.borrow().uff_energy_gradient()
    }
    pub fn uff_energy_gradient_with_params(
        &self,
        params: &ck::UffEvaluationParams,
    ) -> Result<ck::UffEnergyGradient, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.uff_energy_gradient_with_params(params)
        self.inner.borrow().uff_energy_gradient_with_params(params)
    }
    pub fn with_uff_optimized(&self) -> Result<UffOptimizationResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized()
        self.inner
            .borrow()
            .with_uff_optimized()
            .map(|inner| UffOptimizationResult { inner })
    }
    pub fn with_uff_optimized_with_params(
        &self,
        params: &ck::UffOptimizationParams,
    ) -> Result<UffOptimizationResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_with_params(params)
        self.inner
            .borrow()
            .with_uff_optimized_with_params(params)
            .map(|inner| UffOptimizationResult { inner })
    }
    pub fn with_uff_optimized_confs(
        &self,
    ) -> Result<UffConformerOptimizationResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_confs()
        self.inner
            .borrow()
            .with_uff_optimized_confs()
            .map(|inner| UffConformerOptimizationResult { inner })
    }
    pub fn with_uff_optimized_confs_with_params(
        &self,
        params: &ck::UffConformerOptimizationParams,
    ) -> Result<UffConformerOptimizationResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_confs_with_params(params)
        self.inner
            .borrow()
            .with_uff_optimized_confs_with_params(params)
            .map(|inner| UffConformerOptimizationResult { inner })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    const SDF: &str = r#"ethanol
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7000    0.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    1.2000    0.4000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
"#;
    #[test]
    fn uff_all_six_calls_full_values_detached_results_and_errors() {
        let m = Molecule::from_sdf(SDF).unwrap();
        let core = ck::Molecule::from_sdf(SDF).unwrap();
        for p in [
            ck::UffEvaluationParams::default(),
            ck::UffEvaluationParams {
                vdw_threshold: 5.,
                conformer_id: Some(0),
                ignore_interfragment_interactions: false,
            },
        ] {
            let a = m.uff_energy_gradient_with_params(&p).unwrap();
            let b = core.uff_energy_gradient_with_params(&p).unwrap();
            assert_eq!(a.energy().to_bits(), b.energy().to_bits());
            assert_eq!(a.gradient(), b.gradient());
            assert_eq!(a.gradient().len(), 9);
        }
        assert_eq!(
            m.uff_energy_gradient().unwrap().gradient(),
            core.uff_energy_gradient().unwrap().gradient()
        );
        let a = m.with_uff_optimized().unwrap();
        let b = core.with_uff_optimized().unwrap();
        assert_eq!(a.energy().to_bits(), b.energy().to_bits());
        assert_eq!(a.status_code(), b.status_code());
        assert_eq!(a.needs_more(), b.needs_more());
        assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        let p = ck::UffOptimizationParams {
            max_iterations: 2,
            ..Default::default()
        };
        let a = m.with_uff_optimized_with_params(&p).unwrap();
        let b = core.with_uff_optimized_with_params(&p).unwrap();
        assert_eq!(a.energy().to_bits(), b.energy().to_bits());
        assert_eq!(a.status_code(), b.status_code());
        let a = m.with_uff_optimized_confs().unwrap();
        let b = core.with_uff_optimized_confs().unwrap();
        assert_eq!(a.conformer_results(), b.conformer_results());
        assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        let p = ck::UffConformerOptimizationParams {
            max_iterations: 2,
            ..Default::default()
        };
        let a = m.with_uff_optimized_confs_with_params(&p).unwrap();
        let b = core.with_uff_optimized_confs_with_params(&p).unwrap();
        assert_eq!(a.conformer_results(), b.conformer_results());
        assert_eq!(a.conformer_results().len(), 1);
        assert_eq!(*m.inner.borrow(), core);
        let p = ck::UffEvaluationParams {
            conformer_id: Some(99),
            ..Default::default()
        };
        assert!(
            matches!(m.uff_energy_gradient_with_params(&p),Err(ck::OperationError::UffOptimization(ref e)) if e.kind()==ck::UffOptimizationErrorKind::MissingConformer{requested:Some(99)})
        );
        assert!(
            matches!(Molecule::from_smiles("CCO").unwrap().with_uff_optimized(),Err(ck::OperationError::UffOptimization(ref e)) if e.kind()==ck::UffOptimizationErrorKind::MissingConformer{requested:None})
        );
    }
}

#[cfg(test)]
mod multiple_conformer_tests {
    use super::*;
    #[test]
    fn uff_multiple_conformer_order_selectors_and_native_workers() {
        let chemistry = ck::Molecule::from_smiles("CC.CC").unwrap();
        let positions = vec![[0., 0., 0.], [1.9, 0.2, 0.], [5., 1., 0.], [6.7, 1.1, 0.3]];
        let core = ck::Molecule::from_parts(
            chemistry.topology().clone(),
            ck::CoordinateBlock {
                conformers_2d: vec![ck::Conformer2D::new(3, vec![[9., -0.]; 4])],
                conformers_3d: vec![
                    ck::Conformer3D::new(7, positions.clone(), true),
                    ck::Conformer3D::new(3, positions, true),
                ],
                ..Default::default()
            },
            chemistry.properties().clone(),
        )
        .unwrap()
        .with_assigned_valence()
        .unwrap();
        let m = Molecule {
            inner: RefCell::new(core.clone()),
        };
        for threads in [1, 2] {
            let p = ck::UffConformerOptimizationParams {
                num_threads: threads,
                max_iterations: 2,
                ignore_interfragment_interactions: false,
                ..Default::default()
            };
            let a = m.with_uff_optimized_confs_with_params(&p).unwrap();
            let b = core.with_uff_optimized_confs_with_params(&p).unwrap();
            assert_eq!(a.conformer_results(), b.conformer_results());
            assert_eq!(
                a.conformer_results()
                    .iter()
                    .map(|r| r.conformer_id())
                    .collect::<Vec<_>>(),
                [7, 3]
            );
            assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        }
        for id in [None, Some(7), Some(3)] {
            let p = ck::UffOptimizationParams {
                max_iterations: 2,
                conformer_id: id,
                ..Default::default()
            };
            let a = m.with_uff_optimized_with_params(&p).unwrap();
            let b = core.with_uff_optimized_with_params(&p).unwrap();
            assert_eq!(a.energy().to_bits(), b.energy().to_bits());
            assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        }
        assert_eq!(*m.inner.borrow(), core);
    }
}
