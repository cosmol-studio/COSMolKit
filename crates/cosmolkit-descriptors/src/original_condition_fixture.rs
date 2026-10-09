//! Test preparation through the existing parser, sanitize, valence and ring owners.
use crate::{ChiInput, DescriptorInput};
pub(crate) struct Fixture {
    pub(crate) topology: cosmolkit_model::TopologyBlock,
    pub(crate) coordinates: cosmolkit_model::CoordinateBlock,
    pub(crate) properties: cosmolkit_model::MoleculeProperties,
    pub(crate) valence: cosmolkit_core::ValenceAssignment,
    pub(crate) rings: cosmolkit_core::RingInfo,
}
impl Fixture {
    pub(crate) fn from_smiles(smiles: &str) -> Self {
        let parsed = cosmolkit_smiles::parse_smiles(
            smiles,
            &cosmolkit_smiles::SmilesParseParams {
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let removed = cosmolkit_core::remove_hydrogens_with_params(
            parsed.topology,
            parsed.coordinates,
            parsed.properties,
            &cosmolkit_core::RemoveHsParams {
                update_explicit_count: true,
                sanitize: true,
                ..Default::default()
            },
        )
        .unwrap();
        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &removed.topology,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let rings =
            cosmolkit_core::symmetrized_sssr(&removed.topology, &Default::default()).unwrap();
        Self {
            topology: removed.topology,
            coordinates: removed.coordinates,
            properties: removed.properties,
            valence,
            rings,
        }
    }
    pub(crate) fn input(&self) -> DescriptorInput<'_> {
        DescriptorInput::new(
            &self.topology,
            &self.coordinates,
            &self.properties,
            &self.valence,
            &self.rings,
        )
    }
    pub(crate) fn chi(&self) -> ChiInput<'_> {
        ChiInput::new(&self.topology, &self.valence)
    }
}
pub(crate) fn assert_f64_bits(actual: f64, expected: f64, label: &str) {
    assert_eq!(actual.to_bits(), expected.to_bits(), "{label}");
}
pub(crate) fn assert_slice_bits(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(
        actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        "{label}"
    );
}
