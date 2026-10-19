//! Native access to complete-call fingerprint configuration.
//!
//! Composition is owned storage, not a required public nesting hierarchy.
//! Algorithms continue borrowing the canonical generator record directly.

macro_rules! generator_fields {
    ($name:ident, $options:ty, {$($field:ident, $setter:ident: $ty:ty),+ $(,)?}) => {
        impl AsRef<$options> for $name {
            fn as_ref(&self) -> &$options {
                &self.generator
            }
        }
        impl $name {
            $(
                #[must_use]
                pub fn $field(&self) -> $ty {
                    self.generator.$field
                }
                pub fn $setter(&mut self, value: $ty) {
                    self.generator.$field = value;
                }
            )+

            /// Borrow the actual count thresholds owned by this configuration.
            #[must_use]
            pub fn count_bounds(&self) -> &[u32] {
                &self.generator.count_bounds
            }

            /// Edit the actual owned count thresholds without replacing a copy.
            pub fn count_bounds_mut(&mut self) -> &mut Vec<u32> {
                &mut self.generator.count_bounds
            }

            /// Receive independent count thresholds by value.
            pub fn set_count_bounds(&mut self, value: Vec<u32>) {
                self.generator.count_bounds = value;
            }
        }
    };
}
pub(crate) use generator_fields;

#[cfg(test)]
mod tests {
    use crate::{
        AtomPairFingerprintParams, AtomPairGeneratorParams, MorganCallParams,
        MorganFingerprintParams, MorganGeneratorParams, TopologicalTorsionCallParams,
        TopologicalTorsionFingerprintParams, TopologicalTorsionGeneratorParams,
    };

    #[test]
    fn complete_call_configuration_owns_inputs_and_edits_the_native_record() {
        let mut source = MorganGeneratorParams::default();
        source.radius = 2;
        let mut params = MorganFingerprintParams::from(source.clone());
        source.radius = 7;
        source.count_bounds.push(16);
        assert_eq!(params.radius(), 2);
        assert_eq!(params.count_bounds(), &[1, 2, 4, 8]);
        params.set_radius(0);
        params.count_bounds_mut().push(32);
        assert_eq!(params.generator.radius, 0);
        assert_eq!(params.generator.count_bounds, [1, 2, 4, 8, 32]);
        assert_eq!(source.radius, 7);
        assert_eq!(source.count_bounds, [1, 2, 4, 8, 16]);

        let mut pair = AtomPairFingerprintParams::from(AtomPairGeneratorParams::default());
        pair.set_min_distance(2);
        pair.set_include_chirality(true);
        pair.set_count_bounds(vec![]);
        assert_eq!(pair.generator.min_distance, 2);
        assert!(pair.generator.include_chirality);
        assert!(pair.generator.count_bounds.is_empty());

        let mut torsion =
            TopologicalTorsionFingerprintParams::from(TopologicalTorsionGeneratorParams::default());
        torsion.set_torsion_atom_count(3);
        torsion.set_count_simulation(false);
        torsion.count_bounds_mut().clear();
        assert_eq!(torsion.generator.torsion_atom_count, 3);
        assert!(!torsion.generator.count_simulation);
        assert!(torsion.generator.count_bounds.is_empty());
    }

    #[test]
    fn native_composition_moves_buffers_and_preserves_absent_and_empty_selections() {
        let options = MorganGeneratorParams::default();
        let count_bounds = options.count_bounds.as_ptr();
        let call = MorganCallParams::new(Some(vec![]), None, Some(vec![11]), None, 7);
        let atom_invariants = call.custom_atom_invariants.as_ref().unwrap().as_ptr();
        let params = MorganFingerprintParams::from((options, call));
        assert_eq!(params.count_bounds().as_ptr(), count_bounds);
        assert_eq!(params.from_atoms, Some(vec![]));
        assert_eq!(params.ignore_atoms, None);
        assert_eq!(params.conformer_id, 7);
        assert_eq!(
            params.custom_atom_invariants.as_ref().unwrap().as_ptr(),
            atom_invariants
        );

        let call = TopologicalTorsionCallParams {
            ignore_atoms: Some(vec![]),
            use_legacy_stereo_perception: false,
            conformer_id: 9,
            ..Default::default()
        };
        let params = TopologicalTorsionFingerprintParams::from((
            TopologicalTorsionGeneratorParams::default(),
            call,
        ));
        assert_eq!(params.from_atoms, None);
        assert_eq!(params.ignore_atoms, Some(vec![]));
        assert_eq!(params.conformer_id, 9);
        assert!(!params.use_legacy_stereo_perception);
    }
}
