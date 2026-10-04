//! Shared-receiver drawing queries over detached preparation and rendering.

use crate::{DrawingError, Molecule};
use cosmolkit_depict::{DepictOptions, DrawingInput};

#[cfg(all(test, feature = "cap-smiles", feature = "cap-rings"))]
#[path = "drawing_state_probe.rs"]
mod drawing_state_probe;

#[cfg(all(test, feature = "cap-smiles", feature = "cap-rings"))]
#[path = "drawing_svg_boundary_probe.rs"]
mod drawing_svg_boundary_probe;

#[cfg(all(test, feature = "cap-smiles"))]
mod storage_tests {
    use super::*;
    use crate::{Conformer2D, Conformer3D, CoordinateBlock, DerivedState};
    use std::sync::Arc;

    fn coordinate_bits(coordinates: &CoordinateBlock) -> Vec<u64> {
        coordinates
            .conformers_2d
            .iter()
            .flat_map(|c| c.coordinates().iter().flatten().map(|v| v.to_bits()))
            .chain(
                coordinates
                    .conformers_3d
                    .iter()
                    .flat_map(|c| c.coordinates().iter().flatten().map(|v| v.to_bits())),
            )
            .collect()
    }

    #[test]
    fn drawing_storage_eight_queries_preserve_all_four_live_blocks() {
        let parsed = Molecule::from_smiles("CCO").unwrap();
        let coordinates = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(17, vec![[0.0, -0.0], [1.5, 0.0], [2.25, 1.25]])
                    .with_prop("source", "supplied-2d"),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    29,
                    vec![[0.0, 0.0, -0.0], [1.0, 0.5, 0.25], [2.0, 1.0, 1.25]],
                    true,
                )
                .with_prop("source", "independent-3d"),
            ],
            ..CoordinateBlock::default()
        };
        // Test construction retains the parsed, validated cache; no drawing call
        // receives installation authority or writes through these Arc handles.
        let source = Molecule::from_runtime_parts(
            parsed.topology_arc_runtime(),
            Arc::new(coordinates),
            parsed.properties_arc_runtime(),
            parsed.derived_cache_arc_runtime(),
        )
        .unwrap();
        let peer = source.clone();
        let topology = source.topology_arc_runtime();
        let coordinates = source.coordinates_arc_runtime();
        let properties = source.properties_arc_runtime();
        let cache = source.derived_cache_arc_runtime();
        let topology_value = topology.as_ref().clone();
        let coordinate_value = coordinates.as_ref().clone();
        let property_value = properties.as_ref().clone();
        let cache_value = cache.as_ref().clone();
        let bits = coordinate_bits(&coordinates);
        let validity = cache.valid_states();
        assert!(validity.contains(DerivedState::VALENCE));
        let valence = cache.valence_assignment().unwrap() as *const _;
        let valence_value = cache.valence_assignment().unwrap().clone();
        let checkpoint = || {
            for original in [&source, &peer] {
                assert!(Arc::ptr_eq(&topology, &original.topology_arc_runtime()));
                assert!(Arc::ptr_eq(
                    &coordinates,
                    &original.coordinates_arc_runtime()
                ));
                assert!(Arc::ptr_eq(&properties, &original.properties_arc_runtime()));
                assert!(Arc::ptr_eq(&cache, &original.derived_cache_arc_runtime()));
                assert_eq!(original.topology(), &topology_value);
                assert_eq!(original.coordinate_block_runtime(), &coordinate_value);
                assert_eq!(coordinate_bits(original.coordinate_block_runtime()), bits);
                assert_eq!(original.properties(), &property_value);
                assert_eq!(original.derived_cache_runtime(), &cache_value);
                assert_eq!(original.derived_cache_runtime().valid_states(), validity);
                assert_eq!(
                    original
                        .derived_cache_runtime()
                        .valence_assignment()
                        .unwrap() as *const _,
                    valence
                );
                assert_eq!(
                    original.derived_cache_runtime().valence_assignment(),
                    Some(&valence_value)
                );
            }
        };
        let mut calls = 0;
        for _ in 0..2 {
            for molecule in [&source, &peer] {
                for png in [false, true] {
                    if png {
                        checkpoint();
                        let output = molecule.to_png(300, 300);
                        checkpoint();
                        calls += 1;
                        let output = output.unwrap();
                        assert_eq!(&output[..8], b"\x89PNG\r\n\x1a\n");
                    } else {
                        checkpoint();
                        let output = molecule.to_svg(300, 300);
                        checkpoint();
                        calls += 1;
                        assert!(output.unwrap().contains("</svg>"));
                    }
                }
            }
        }
        assert_eq!(calls, 8);
    }
}

impl Molecule {
    fn drawing_input(&self) -> DrawingInput<'_> {
        #[cfg(any(
            feature = "cap-valence",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors"
        ))]
        let valence = self.derived_cache_runtime().valence_assignment();
        #[cfg(not(any(
            feature = "cap-valence",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors"
        )))]
        let valence = None;
        #[cfg(feature = "cap-rings")]
        let rings = self.derived_cache_runtime().valid_ring_info();
        #[cfg(not(feature = "cap-rings"))]
        let rings = None;
        DrawingInput {
            topology: self.topology(),
            coordinates: self.coordinate_block_runtime(),
            properties: self.properties(),
            valence,
            rings,
        }
    }

    /// Render Experimental SVG using the first stored 2D layout, or generate
    /// a detached layout when absent. Stored coordinates and caches are preserved.
    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, DrawingError> {
        cosmolkit_depict::render_svg(self.drawing_input(), &DepictOptions { width, height })
    }

    /// Rasterize the same Experimental SVG with the embedded Noto Sans font.
    /// Preparation and rendering never install changes in this molecule.
    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, DrawingError> {
        cosmolkit_depict::render_png(self.drawing_input(), &DepictOptions { width, height })
    }
}
