//! Thin MMFF evaluation and detached optimization projections.
use crate::Molecule;
use cosmolkit as ck;
use std::cell::RefCell;
pub struct MmffOptimizeMoleculeResult {
    inner: ck::MmffOptimizeMoleculeResult,
}
impl MmffOptimizeMoleculeResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule.clone(),
        Molecule {
            inner: RefCell::new(self.inner.molecule.clone()),
        }
    }
    pub fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    pub fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
}
pub struct MmffOptimizeMoleculeConfsResult {
    inner: ck::MmffOptimizeMoleculeConfsResult,
}
impl MmffOptimizeMoleculeConfsResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule.clone(),
        Molecule {
            inner: RefCell::new(self.inner.molecule.clone()),
        }
    }
    pub fn conformer_results(&self) -> &[ck::MmffOptimizeMoleculeConfResult] {
        &self.inner.conformer_results
    }
}
impl Molecule {
    pub fn mmff_energy_gradient(
        &self,
    ) -> Result<Option<ck::MmffEnergyGradient>, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.mmff_energy_gradient()
        self.inner.borrow().mmff_energy_gradient()
    }
    pub fn mmff_energy_gradient_with_params(
        &self,
        params: &ck::MmffEvaluationParams,
    ) -> Result<Option<ck::MmffEnergyGradient>, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.mmff_energy_gradient_with_params(params)
        self.inner.borrow().mmff_energy_gradient_with_params(params)
    }
    pub fn with_mmff_optimized(&self) -> Result<MmffOptimizeMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized()
        self.inner
            .borrow()
            .with_mmff_optimized()
            .map(|inner| MmffOptimizeMoleculeResult { inner })
    }
    pub fn with_mmff_optimized_with_params(
        &self,
        params: &ck::MmffOptimizationParams,
    ) -> Result<MmffOptimizeMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_with_params(params)
        self.inner
            .borrow()
            .with_mmff_optimized_with_params(params)
            .map(|inner| MmffOptimizeMoleculeResult { inner })
    }
    pub fn with_mmff_optimized_conformers(
        &self,
    ) -> Result<MmffOptimizeMoleculeConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_confs()
        self.inner
            .borrow()
            .with_mmff_optimized_conformers()
            .map(|inner| MmffOptimizeMoleculeConfsResult { inner })
    }
    pub fn with_mmff_optimized_conformers_with_params(
        &self,
        params: &ck::MmffConformerOptimizationParams,
    ) -> Result<MmffOptimizeMoleculeConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_confs_with_params(params)
        self.inner
            .borrow()
            .with_mmff_optimized_conformers_with_params(params)
            .map(|inner| MmffOptimizeMoleculeConfsResult { inner })
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
    fn mmff_all_six_calls_variants_values_and_source_invalid_type_outcomes() {
        let m = Molecule::from_sdf(SDF).unwrap();
        let core = ck::Molecule::from_sdf(SDF).unwrap();
        assert_eq!(
            m.mmff_energy_gradient().unwrap(),
            core.mmff_energy_gradient().unwrap()
        );
        for variant in ["MMFF94", "MMFF94s"] {
            let p = ck::MmffEvaluationParams {
                mmff_variant: variant.into(),
                non_bonded_threshold: 5.,
                conformer_id: Some(0),
                ignore_interfragment_interactions: false,
            };
            let a = m.mmff_energy_gradient_with_params(&p).unwrap().unwrap();
            let b = core.mmff_energy_gradient_with_params(&p).unwrap().unwrap();
            assert_eq!(a.energy().to_bits(), b.energy().to_bits());
            assert_eq!(a.gradient(), b.gradient());
            assert_eq!(a.gradient().len(), 9);
            let p = ck::MmffOptimizationParams {
                mmff_variant: variant.into(),
                max_iterations: 2,
                ..Default::default()
            };
            let a = m.with_mmff_optimized_with_params(&p).unwrap();
            let b = core.with_mmff_optimized_with_params(&p).unwrap();
            assert_eq!(a.status_code(), b.status_code());
            assert_eq!(a.needs_more(), b.needs_more());
            assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
            let p = ck::MmffConformerOptimizationParams {
                mmff_variant: variant.into(),
                max_iterations: 2,
                ..Default::default()
            };
            let a = m.with_mmff_optimized_conformers_with_params(&p).unwrap();
            let b = core.with_mmff_optimized_conformers_with_params(&p).unwrap();
            assert_eq!(a.conformer_results(), b.conformer_results());
            assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        }
        let a = m.with_mmff_optimized().unwrap();
        let b = core.with_mmff_optimized().unwrap();
        assert_eq!(a.status_code(), b.status_code());
        assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        let a = m.with_mmff_optimized_conformers().unwrap();
        let b = core.with_mmff_optimized_conformers().unwrap();
        assert_eq!(a.conformer_results(), b.conformer_results());
        assert_eq!(*m.inner.borrow(), core);
        let dummy = Molecule::from_smiles("*").unwrap();
        assert_eq!(dummy.mmff_energy_gradient().unwrap(), None);
        assert_eq!(dummy.with_mmff_optimized().unwrap().status_code(), -1);
        assert!(
            dummy
                .with_mmff_optimized_conformers()
                .unwrap()
                .conformer_results()
                .is_empty()
        );
        assert!(matches!(
            Molecule::from_smiles("CCO").unwrap().mmff_energy_gradient(),
            Err(ck::OperationError::MmffOptimization(_))
        ));
    }
    #[test]
    fn mmff_multiple_conformer_order_and_native_workers() {
        let chemistry = ck::Molecule::from_smiles("CCO").unwrap();
        let positions = vec![[0., 0., 0.], [1.7, 0.2, 0.], [2.5, 1.2, 0.4]];
        let core = ck::Molecule::from_parts(
            chemistry.topology().clone(),
            ck::CoordinateBlock {
                conformers_3d: vec![
                    ck::Conformer3D::new(7, positions.clone(), true),
                    ck::Conformer3D::new(3, positions, true),
                ],
                ..Default::default()
            },
            chemistry.properties().clone(),
        )
        .unwrap();
        let m = Molecule {
            inner: RefCell::new(core.clone()),
        };
        for threads in [1, 2] {
            let p = ck::MmffConformerOptimizationParams {
                num_threads: threads,
                max_iterations: 2,
                ..Default::default()
            };
            let a = m.with_mmff_optimized_conformers_with_params(&p).unwrap();
            let b = core.with_mmff_optimized_conformers_with_params(&p).unwrap();
            assert_eq!(a.conformer_results(), b.conformer_results());
            assert_eq!(a.conformer_results().len(), 2);
            assert_eq!(*a.molecule().inner.borrow(), *b.molecule());
        }
        for id in [None, Some(7), Some(3)] {
            let p = ck::MmffEvaluationParams {
                conformer_id: id,
                ..Default::default()
            };
            assert_eq!(
                m.mmff_energy_gradient_with_params(&p).unwrap(),
                core.mmff_energy_gradient_with_params(&p).unwrap()
            );
        }
        assert_eq!(*m.inner.borrow(), core);
    }
}
