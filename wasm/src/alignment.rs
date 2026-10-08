//! Thin alignment methods; only the canonical public facade performs chemistry.
use super::Molecule;
use cosmolkit as ck;

// Host-wrapper policy from the pinned Python projection. This is conversion,
// not alignment logic: empty optional lists are represented as absent pointers.
fn projected_alignment(
    params: &ck::AlignmentParameters,
    atom_count: usize,
) -> Result<ck::AlignmentParameters, ck::AlignmentError> {
    // RDKit❗✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit❗✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit❗✔️:   unsigned int nDoubles = doubles.size();
    // RDKit❗✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit❗✔️:   doubleVec = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   if (nDoubles > 0) {
    // RDKit❗✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit❗✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit❗✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return doubleVec;
    // RDKit❗✔️: }
    // RDKit❗✔️: std::vector<std::pair<int, int>> *translateAtomMap(
    // RDKit❗✔️:     const python::object &atomMap) {
    // RDKit❗✔️:   PySequenceHolder<python::object> pyAtomMap(atomMap);
    // RDKit❗✔️:   std::vector<std::pair<int, int>> *res;
    // RDKit❗✔️:   res = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   unsigned int n = pyAtomMap.size();
    // RDKit❗✔️:   if (n > 0) {
    // RDKit❗✔️:     res = new std::vector<std::pair<int, int>>;
    // RDKit❗✔️:     for (i = 0; i < n; ++i) {
    // RDKit❗✔️:       PySequenceHolder<int> item(pyAtomMap[i]);
    // RDKit❗✔️:       if (item.size() != 2) {
    // RDKit❗✔️:         delete res;
    // RDKit❗✔️:         res = nullptr;
    // RDKit❗✔️:         throw_value_error("Incorrect format for an atomMap");
    // RDKit❗✔️:       }
    // RDKit❗✔️:       res->push_back(std::pair<int, int>(item[0], item[1]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let mut projected = params.clone();
    projected.atom_map = projected.atom_map.filter(|v| !v.is_empty());
    projected.weights = projected.weights.filter(|v| !v.is_empty());
    let expected = projected.atom_map.as_ref().map_or(atom_count, Vec::len);
    check_weight_count(projected.weights.as_deref(), expected)?;
    Ok(projected)
}
fn check_weight_count(weights: Option<&[f64]>, expected: usize) -> Result<(), ck::AlignmentError> {
    // RDKit❗✔️:   if (wtsVec) {
    // RDKit❗✔️:     if (wtsVec->size() != nAtms) {
    // RDKit❗✔️:       throw_value_error("Incorrect number of weights specified");
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    if let Some(weights) = weights {
        if weights.len() != expected {
            return Err(ck::AlignmentError::WeightCountMismatch {
                map_len: expected,
                weight_len: weights.len(),
            });
        }
    }
    Ok(())
}
fn projected_best(
    params: &ck::BestAlignmentParameters,
) -> Result<ck::BestAlignmentParameters, ck::AlignmentError> {
    // RDKit❗✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit❗✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit❗✔️:   unsigned int nDoubles = doubles.size();
    // RDKit❗✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit❗✔️:   doubleVec = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   if (nDoubles > 0) {
    // RDKit❗✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit❗✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit❗✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return doubleVec;
    // RDKit❗✔️: }
    let mut projected = params.clone();
    projected.weights = projected.weights.filter(|v| !v.is_empty());
    if let Some(first) = projected.atom_maps.first() {
        check_weight_count(projected.weights.as_deref(), first.len())?;
    }
    Ok(projected)
}
fn projected_all(
    params: &ck::AllConformerRmsdParameters,
) -> Result<ck::AllConformerRmsdParameters, ck::AlignmentError> {
    // RDKit❗✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit❗✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit❗✔️:   unsigned int nDoubles = doubles.size();
    // RDKit❗✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit❗✔️:   doubleVec = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   if (nDoubles > 0) {
    // RDKit❗✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit❗✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit❗✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return doubleVec;
    // RDKit❗✔️: }
    let mut projected = params.clone();
    projected.weights = projected.weights.filter(|v| !v.is_empty());
    if let Some(first) = projected.atom_maps.first() {
        check_weight_count(projected.weights.as_deref(), first.len())?;
    }
    Ok(projected)
}
fn projected_coordinate(params: &ck::CoordinateRmsdParameters) -> ck::CoordinateRmsdParameters {
    // RDKit❗✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit❗✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit❗✔️:   unsigned int nDoubles = doubles.size();
    // RDKit❗✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit❗✔️:   doubleVec = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   if (nDoubles > 0) {
    // RDKit❗✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit❗✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit❗✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return doubleVec;
    // RDKit❗✔️: }
    let mut projected = params.clone();
    projected.weights = projected.weights.filter(|v| !v.is_empty());
    projected
}
fn projected_conformer(
    params: &ck::ConformerAlignmentParameters,
) -> ck::ConformerAlignmentParameters {
    // RDKit❗✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit❗✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit❗✔️:   unsigned int nDoubles = doubles.size();
    // RDKit❗✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit❗✔️:   doubleVec = nullptr;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   if (nDoubles > 0) {
    // RDKit❗✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit❗✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit❗✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return doubleVec;
    // RDKit❗✔️: }
    let mut projected = params.clone();
    projected.weights = projected.weights.filter(|v| !v.is_empty());
    projected
}

impl Molecule {
    pub fn alignment_transform_to(
        &self,
        reference: &Molecule,
    ) -> Result<ck::AlignmentResult, ck::AlignmentError> {
        self.inner
            .borrow()
            .alignment_transform_to(&reference.inner.borrow())
    }
    pub fn alignment_transform_to_with_params(
        &self,
        reference: &Molecule,
        params: &ck::AlignmentParameters,
    ) -> Result<ck::AlignmentResult, ck::AlignmentError> {
        let params = projected_alignment(params, self.inner.borrow().num_atoms())?;
        self.inner
            .borrow()
            .alignment_transform_to_with_params(&reference.inner.borrow(), &params)
    }
    pub fn best_alignment_to(
        &self,
        reference: &Molecule,
    ) -> Result<ck::AlignmentResult, ck::AlignmentError> {
        self.inner
            .borrow()
            .best_alignment_to(&reference.inner.borrow())
    }
    pub fn best_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &ck::BestAlignmentParameters,
    ) -> Result<ck::AlignmentResult, ck::AlignmentError> {
        let params = projected_best(params)?;
        self.inner
            .borrow()
            .best_alignment_to_with_params(&reference.inner.borrow(), &params)
    }
    pub fn best_rmsd_to(&self, reference: &Molecule) -> Result<f64, ck::AlignmentError> {
        self.inner.borrow().best_rmsd_to(&reference.inner.borrow())
    }
    pub fn best_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &ck::BestAlignmentParameters,
    ) -> Result<f64, ck::AlignmentError> {
        let params = projected_best(params)?;
        self.inner
            .borrow()
            .best_rmsd_to_with_params(&reference.inner.borrow(), &params)
    }
    pub fn coordinate_rmsd_to(&self, reference: &Molecule) -> Result<f64, ck::AlignmentError> {
        self.inner
            .borrow()
            .coordinate_rmsd_to(&reference.inner.borrow())
    }
    pub fn coordinate_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &ck::CoordinateRmsdParameters,
    ) -> Result<f64, ck::AlignmentError> {
        let params = projected_coordinate(params);
        self.inner
            .borrow()
            .coordinate_rmsd_to_with_params(&reference.inner.borrow(), &params)
    }
    pub fn all_conformer_best_rmsds(&self) -> Result<Vec<ck::ConformerRmsd>, ck::AlignmentError> {
        self.inner.borrow().all_conformer_best_rmsds()
    }
    pub fn all_conformer_best_rmsds_with_params(
        &self,
        params: &ck::AllConformerRmsdParameters,
    ) -> Result<Vec<ck::ConformerRmsd>, ck::AlignmentError> {
        let params = projected_all(params)?;
        self.inner
            .borrow()
            .all_conformer_best_rmsds_with_params(&params)
    }
    pub fn with_alignment_to(
        &self,
        reference: &Molecule,
    ) -> Result<(Molecule, ck::AlignmentResult), ck::OperationError> {
        self.inner
            .borrow()
            .with_alignment_to(&reference.inner.borrow())
            .map(|(inner, report)| {
                (
                    Molecule {
                        inner: inner.into(),
                    },
                    report,
                )
            })
    }
    pub fn with_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &ck::AlignmentParameters,
    ) -> Result<(Molecule, ck::AlignmentResult), ck::OperationError> {
        let params = projected_alignment(params, self.inner.borrow().num_atoms())
            .map_err(ck::OperationError::Alignment)?;
        self.inner
            .borrow()
            .with_alignment_to_with_params(&reference.inner.borrow(), &params)
            .map(|(inner, report)| {
                (
                    Molecule {
                        inner: inner.into(),
                    },
                    report,
                )
            })
    }
    pub fn align_to_(
        &self,
        reference: &Molecule,
    ) -> Result<ck::AlignmentResult, ck::OperationError> {
        let reference = reference.inner.borrow().clone();
        self.inner.borrow_mut().align_to_(&reference)
    }
    pub fn align_to_with_params_(
        &self,
        reference: &Molecule,
        params: &ck::AlignmentParameters,
    ) -> Result<ck::AlignmentResult, ck::OperationError> {
        let params = projected_alignment(params, self.inner.borrow().num_atoms())
            .map_err(ck::OperationError::Alignment)?;
        let reference = reference.inner.borrow().clone();
        self.inner
            .borrow_mut()
            .align_to_with_params_(&reference, &params)
    }
    pub fn with_aligned_conformers(
        &self,
    ) -> Result<(Molecule, ck::ConformerAlignmentReport), ck::OperationError> {
        self.inner
            .borrow()
            .with_aligned_conformers()
            .map(|(inner, report)| {
                (
                    Molecule {
                        inner: inner.into(),
                    },
                    report,
                )
            })
    }
    pub fn with_aligned_conformers_with_params(
        &self,
        params: &ck::ConformerAlignmentParameters,
    ) -> Result<(Molecule, ck::ConformerAlignmentReport), ck::OperationError> {
        let params = projected_conformer(params);
        self.inner
            .borrow()
            .with_aligned_conformers_with_params(&params)
            .map(|(inner, report)| {
                (
                    Molecule {
                        inner: inner.into(),
                    },
                    report,
                )
            })
    }
    pub fn align_conformers_(&self) -> Result<ck::ConformerAlignmentReport, ck::OperationError> {
        self.inner.borrow_mut().align_conformers_()
    }
    pub fn align_conformers_with_params_(
        &self,
        params: &ck::ConformerAlignmentParameters,
    ) -> Result<ck::ConformerAlignmentReport, ck::OperationError> {
        let params = projected_conformer(params);
        self.inner
            .borrow_mut()
            .align_conformers_with_params_(&params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source(shift: f64) -> Molecule {
        let molecule = ck::Molecule::from_smiles("CO").unwrap();
        let mut builder = molecule.to_builder();
        builder
            .add_3d_conformer(vec![[shift, 0., 0.], [shift + 1., 0., 0.]])
            .unwrap();
        builder
            .add_3d_conformer(vec![[shift + 4., 0., 0.], [shift + 5., 0., 0.]])
            .unwrap();
        Molecule {
            inner: builder.build().unwrap().into(),
        }
    }

    #[test]
    fn alignment_projection_preserves_results_value_returns_and_self_aliasing() {
        let probe = source(3.);
        let reference = source(0.);
        let result = probe.alignment_transform_to(&reference).unwrap();
        assert_eq!(result.atom_map().len(), 2);
        assert!(result.rmsd() < 1e-12);
        assert!((result.transform().matrix()[0][3] + 3.).abs() < 1e-12);
        let original = probe.coordinates_3d(0);
        let (aligned, _) = probe.with_alignment_to(&reference).unwrap();
        assert_eq!(probe.coordinates_3d(0), original);
        assert!(aligned.coordinate_rmsd_to(&reference).unwrap() < 1e-12);
        assert_eq!(probe.all_conformer_best_rmsds().unwrap().len(), 1);
        let pair = probe.all_conformer_best_rmsds().unwrap()[0];
        assert_eq!(pair.probe_conformer_id(), 1);
        assert_eq!(pair.reference_conformer_id(), 0);
        assert!(pair.rmsd() < 1e-12);
        let (_, report) = probe.with_aligned_conformers().unwrap();
        assert_eq!(report.rmsds().len(), 1);
        assert!(probe.align_to_(&probe).unwrap().rmsd() < 1e-12);
        assert_eq!(probe.coordinates_3d(0), original);
        assert_eq!(probe.align_conformers_().unwrap().rmsds().len(), 1);
    }

    #[test]
    fn alignment_projection_keeps_empty_pointer_policy_and_typed_failures() {
        let probe = source(3.);
        let reference = source(0.);
        let params = ck::AlignmentParameters {
            atom_map: Some(vec![]),
            weights: Some(vec![]),
            ..Default::default()
        };
        assert!(
            probe
                .alignment_transform_to_with_params(&reference, &params)
                .unwrap()
                .rmsd()
                < 1e-12
        );
        assert_eq!(params.atom_map, Some(vec![]));
        assert_eq!(params.weights, Some(vec![]));
        let params = ck::AlignmentParameters {
            weights: Some(vec![1.]),
            ..Default::default()
        };
        assert!(matches!(
            probe.alignment_transform_to_with_params(&reference, &params),
            Err(ck::AlignmentError::WeightCountMismatch {
                map_len: 2,
                weight_len: 1
            })
        ));
        let before = probe.coordinates_3d(0);
        assert!(matches!(
            probe.align_to_with_params_(&reference, &params),
            Err(ck::OperationError::Alignment(
                ck::AlignmentError::WeightCountMismatch {
                    map_len: 2,
                    weight_len: 1
                }
            ))
        ));
        assert_eq!(probe.coordinates_3d(0), before);
        let params = ck::AlignmentParameters {
            probe_conformer_id: 77,
            ..Default::default()
        };
        assert!(matches!(
            probe.alignment_transform_to_with_params(&reference, &params),
            Err(ck::AlignmentError::ConformerNotFound { id: 77 })
        ));
    }
}
