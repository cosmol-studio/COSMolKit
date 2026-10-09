//! Thin typed user-property operations over the canonical detached Atom store.
#[cfg(feature = "cap-transforms")]
use super::PreservationProof;
use crate::{AtomId, Molecule, OperationError, PropertyValue};
#[cfg(feature = "cap-transforms")]
use crate::{DerivedState, TopologyEditKind};
#[cfg(feature = "cap-transforms")]
use cosmolkit_macros::mol_op_body;

impl Molecule {
    /// Read a property without string conversion. Missing keys return `None`;
    /// invalid atom IDs return a structured error. This never recalculates state.
    pub fn atom_property(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<&PropertyValue>, OperationError> {
        let row = self.atom(atom).ok_or(OperationError::AtomPropertyIndex {
            atom,
            atom_count: self.num_atoms(),
        })?;
        Ok(row.prop(key))
    }
}

#[cfg(feature = "cap-transforms")]
#[mol_op_body(with_atom_property, parts)]
pub(crate) fn with_atom_property_impl(
    atom: AtomId,
    key: &str,
    value: &PropertyValue,
) -> Result<(), OperationError> {
    if !cosmolkit_model::is_user_atom_property(key.as_bytes()) {
        return Err(OperationError::ReservedAtomPropertyKey {
            key: key.to_owned(),
        });
    }
    let mut topology = parts.checkout_topology()?;
    // Always restore the owned block before propagating any fallible result.
    // MODEL validates the key before writing. No coordinate/property/cache clone
    // is used to implement rollback; the source value remains protected by COW.
    let result = (|| {
        let count = topology.atoms.len();
        let row =
            topology
                .atoms
                .get_mut(atom.index())
                .ok_or(OperationError::AtomPropertyIndex {
                    atom,
                    atom_count: count,
                })?;
        if row
            .is_prop_computed(key)
            .map_err(cosmolkit_model::AtomPropertyError::ComputedListKind)
            .map_err(OperationError::AtomProperty)?
        {
            return Err(OperationError::ReservedAtomPropertyKey {
                key: key.to_owned(),
            });
        }
        row.set_prop(key, value)
            .map_err(OperationError::AtomProperty)
    })();
    parts.install_topology(topology)?;
    result?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::DRAWING.union(DerivedState::FINGERPRINT))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES),
        PreservationProof::AtomUserProperties,
    )?;
    parts.apply_cip_policy()
}

#[cfg(all(test, feature = "cap-smiles", feature = "cap-transforms"))]
mod tests {
    use super::*;
    #[test]
    fn typed_atom_property_value_inplace_cow_and_identity() {
        let source = Molecule::from_smiles("CCO").unwrap();
        let peer = source.clone();
        let atom = AtomId::new(0);
        for value in [
            PropertyValue::Int(42),
            PropertyValue::UInt(u32::MAX),
            PropertyValue::Bool(true),
            PropertyValue::Double(-0.0),
            PropertyValue::from("note"),
            PropertyValue::IntVector(vec![1, 2]),
            PropertyValue::StringVector(vec!["a".into(), "b".into()]),
        ] {
            let mut tagged = source
                .with_atom_property(atom, "tracking_id", &value)
                .unwrap();
            assert_eq!(
                tagged.atom_property(atom, "tracking_id").unwrap(),
                Some(&value)
            );
            assert_eq!(tagged.to_smiles().unwrap(), source.to_smiles().unwrap());
            assert!(std::ptr::eq(
                tagged.coordinate_block_runtime(),
                source.coordinate_block_runtime()
            ));
            assert!(std::ptr::eq(tagged.properties(), source.properties()));
            assert!(source.atom_property(atom, "tracking_id").unwrap().is_none());
            assert!(peer.atom_property(atom, "tracking_id").unwrap().is_none());
            let copied = tagged.clone();
            tagged
                .set_atom_property_(atom, "tracking_id", &PropertyValue::Int(7))
                .unwrap();
            assert_eq!(
                copied.atom_property(atom, "tracking_id").unwrap(),
                Some(&value)
            );
            assert_eq!(
                tagged.atom_property(atom, "tracking_id").unwrap(),
                Some(&PropertyValue::Int(7))
            );
        }
    }
    #[test]
    fn typed_atom_property_failures_preserve_receiver_and_peer() {
        let mut mol = Molecule::from_smiles("CCO").unwrap();
        let peer = mol.clone();
        for key in [
            "",
            "_CIPCode",
            "__computedProps",
            "molAtomMapNumber",
            "react_idx",
        ] {
            assert!(matches!(
                mol.set_atom_property_(AtomId::new(0), key, &42.into()),
                Err(OperationError::ReservedAtomPropertyKey { .. })
            ));
        }
        assert!(matches!(
            mol.set_atom_property_(AtomId::new(99), "tracking_id", &42.into()),
            Err(OperationError::AtomPropertyIndex { .. })
        ));
        assert!(matches!(
            mol.atom_property(AtomId::new(99), "tracking_id"),
            Err(OperationError::AtomPropertyIndex { .. })
        ));
        assert_eq!(mol.topology(), peer.topology());
        assert_eq!(mol.properties(), peer.properties());
    }
    #[cfg(feature = "cap-transforms")]
    #[test]
    fn typed_atom_property_follows_fragment_atoms() {
        let source = Molecule::from_smiles("CC.O")
            .unwrap()
            .with_atom_property(AtomId::new(2), "tracking_id", &84.into())
            .unwrap();
        let fragments = source.fragments().unwrap();
        assert_eq!(fragments.len(), 2);
        assert!(
            fragments[0]
                .atom_property(AtomId::new(0), "tracking_id")
                .unwrap()
                .is_none()
        );
        assert_eq!(
            fragments[1]
                .atom_property(AtomId::new(0), "tracking_id")
                .unwrap(),
            Some(&PropertyValue::Int(84))
        );
        assert_eq!(
            source.atom_property(AtomId::new(2), "tracking_id").unwrap(),
            Some(&PropertyValue::Int(84))
        );
    }

    #[test]
    fn typed_atom_property_preserves_and_rejects_computed_user_keys() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_computed_prop("derived_score", 7)
                .unwrap(),
        );
        let mut mol = builder.build().unwrap();
        let peer = mol.clone();
        let atom = AtomId::new(0);
        assert!(matches!(
            mol.set_atom_property_(atom, "derived_score", &42.into()),
            Err(OperationError::ReservedAtomPropertyKey { .. })
        ));
        assert_eq!(mol.topology(), peer.topology());
        let tagged = mol
            .with_atom_property(atom, "tracking_id", &42.into())
            .unwrap();
        assert_eq!(
            tagged.atom_property(atom, "derived_score").unwrap(),
            Some(&7.into())
        );
        assert!(
            tagged
                .atom(atom)
                .unwrap()
                .is_prop_computed("derived_score")
                .unwrap()
        );
        assert!(mol.atom_property(atom, "tracking_id").unwrap().is_none());
    }
}
