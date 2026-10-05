//! Thin canonical alignment APIs; all algorithms live in the detached owner.
use crate::{
    AlignmentError, AlignmentParameters, AlignmentResult, AllConformerRmsdParameters,
    BestAlignmentParameters, ConformerRmsd, CoordinateRmsdParameters, Molecule,
};
use cosmolkit_alignment::AlignmentInput;

pub(crate) fn input<'a>(
    topology: &'a crate::TopologyBlock,
    coordinates: &'a crate::CoordinateBlock,
    cache: &'a crate::molecule::DerivedCacheBlock,
) -> AlignmentInput<'a> {
    #[cfg(any(
        feature = "cap-valence",
        feature = "cap-fingerprints",
        feature = "cap-hydrogens",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-descriptors",
        feature = "cap-forcefields",
        feature = "cap-tautomer",
        feature = "cap-hashing",
        feature = "cap-serialization"
    ))]
    let valence = cache.valence_assignment();
    #[cfg(not(any(
        feature = "cap-valence",
        feature = "cap-fingerprints",
        feature = "cap-hydrogens",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-descriptors",
        feature = "cap-forcefields",
        feature = "cap-tautomer",
        feature = "cap-hashing",
        feature = "cap-serialization"
    )))]
    let valence = None;
    #[cfg(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-forcefields",
        feature = "cap-fingerprints",
        feature = "cap-tautomer",
        feature = "cap-hashing",
        feature = "cap-serialization"
    ))]
    let rings = cache.valid_ring_info();
    #[cfg(not(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-forcefields",
        feature = "cap-fingerprints",
        feature = "cap-tautomer",
        feature = "cap-hashing",
        feature = "cap-serialization"
    )))]
    let rings = None;
    let _ = cache;
    AlignmentInput {
        topology,
        coordinates,
        valence,
        rings,
    }
}
impl Molecule {
    fn alignment_input(&self) -> AlignmentInput<'_> {
        input(
            self.topology(),
            self.coordinate_block_runtime(),
            self.derived_cache_runtime(),
        )
    }
    pub fn alignment_transform_to(
        &self,
        reference: &Self,
    ) -> Result<AlignmentResult, AlignmentError> {
        self.alignment_transform_to_with_params(reference, &AlignmentParameters::default())
    }
    pub fn alignment_transform_to_with_params(
        &self,
        reference: &Self,
        params: &AlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        self.alignment_input()
            .alignment_transform_to(&reference.alignment_input(), params)
    }
    pub fn best_alignment_to(&self, reference: &Self) -> Result<AlignmentResult, AlignmentError> {
        self.best_alignment_to_with_params(reference, &BestAlignmentParameters::default())
    }
    pub fn best_alignment_to_with_params(
        &self,
        reference: &Self,
        params: &BestAlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        self.alignment_input()
            .best_alignment_to(&reference.alignment_input(), params)
    }
    pub fn best_rmsd_to(&self, reference: &Self) -> Result<f64, AlignmentError> {
        self.best_rmsd_to_with_params(reference, &BestAlignmentParameters::default())
    }
    pub fn best_rmsd_to_with_params(
        &self,
        reference: &Self,
        params: &BestAlignmentParameters,
    ) -> Result<f64, AlignmentError> {
        self.alignment_input()
            .best_rmsd_to(&reference.alignment_input(), params)
    }
    pub fn coordinate_rmsd_to(&self, reference: &Self) -> Result<f64, AlignmentError> {
        self.coordinate_rmsd_to_with_params(reference, &CoordinateRmsdParameters::default())
    }
    pub fn coordinate_rmsd_to_with_params(
        &self,
        reference: &Self,
        params: &CoordinateRmsdParameters,
    ) -> Result<f64, AlignmentError> {
        self.alignment_input()
            .coordinate_rmsd_to(&reference.alignment_input(), params)
    }
    pub fn all_conformer_best_rmsds(&self) -> Result<Vec<ConformerRmsd>, AlignmentError> {
        self.all_conformer_best_rmsds_with_params(&AllConformerRmsdParameters::default())
    }
    pub fn all_conformer_best_rmsds_with_params(
        &self,
        params: &AllConformerRmsdParameters,
    ) -> Result<Vec<ConformerRmsd>, AlignmentError> {
        self.alignment_input().all_conformer_best_rmsds(params)
    }
}

pub(crate) fn alignment_transform_candidate(
    topology: &crate::TopologyBlock,
    coordinates: &crate::CoordinateBlock,
    cache: &crate::molecule::DerivedCacheBlock,
    reference: &Molecule,
    params: &AlignmentParameters,
) -> Result<AlignmentResult, AlignmentError> {
    input(topology, coordinates, cache).alignment_transform_to(&reference.alignment_input(), params)
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    use super::*;
    use std::sync::Arc;
    fn source(shift: f64) -> Molecule {
        let m = Molecule::from_smiles("CC").unwrap();
        let mut b = m.to_builder();
        b.add_3d_conformer(vec![[shift, 0., 0.], [shift + 1., 0., 0.]])
            .unwrap();
        b.build().unwrap()
    }
    fn unchanged(a: &Molecule, b: &Molecule) {
        assert_eq!(a, b);
        assert!(Arc::ptr_eq(
            &a.topology_arc_runtime(),
            &b.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &a.coordinates_arc_runtime(),
            &b.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &a.properties_arc_runtime(),
            &b.properties_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &a.derived_cache_arc_runtime(),
            &b.derived_cache_arc_runtime()
        ));
    }
    #[test]
    fn readonly_and_value_alignment_preserve_source_blocks() {
        let m = source(3.);
        let peer = m.clone();
        let r = source(0.);
        let rp = r.clone();
        m.best_alignment_to(&r).unwrap();
        unchanged(&m, &peer);
        let (a, report) = m.with_alignment_to(&r).unwrap();
        assert_eq!(report.rmsd(), 0.);
        unchanged(&m, &peer);
        unchanged(&r, &rp);
        assert!(Arc::ptr_eq(
            &m.topology_arc_runtime(),
            &a.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &m.properties_arc_runtime(),
            &a.properties_arc_runtime()
        ));
        assert_eq!(m.derived_cache_runtime(), a.derived_cache_runtime());
        assert!(!Arc::ptr_eq(
            &m.coordinates_arc_runtime(),
            &a.coordinates_arc_runtime()
        ));
    }
    #[test]
    fn failed_alignment_returns_no_report_and_preserves_every_block() {
        let mut m = source(3.);
        let peer = m.clone();
        let r = source(0.);
        let p = AlignmentParameters {
            atom_map: Some(vec![crate::AlignmentAtomMap::new(99, 0)]),
            ..Default::default()
        };
        assert!(m.align_to_with_params_(&r, &p).is_err());
        unchanged(&m, &peer);
        assert!(m.with_alignment_to_with_params(&r, &p).is_err());
        unchanged(&m, &peer);
    }
}
