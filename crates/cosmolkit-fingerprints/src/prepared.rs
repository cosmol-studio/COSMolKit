//! Source-ordered preparation of the temporary Morgan environment molecule.
//!
//! Source attribution: pinned RDKit `GraphMol/Fingerprints/
//! FingerprintGenerator.cpp`, `GraphMol/ROMol.{cpp,h}`, and the
//! stereochemistry assignment source declared in `THIRD_PARTY_NOTICES.md`.

use std::borrow::Cow;

use cosmolkit_core::{RingInfo, ValenceAssignment, assign_legacy_stereochemistry_for_depiction};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use crate::MorganError;

/// Borrowed final molecule state and caller-prepared assignments for one
/// detached Morgan fingerprint call.
#[derive(Debug, Clone, Copy)]
pub struct MorganPreparedInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub valence: &'a ValenceAssignment,
    pub rings: &'a RingInfo,
}

/// Environment view matching RDKit's conditional `lmol` selection.
///
/// The original input remains owned by the caller and is still the source for
/// atom/bond invariant generation. Only the environment sees this view.
#[derive(Debug)]
pub(crate) struct PreparedMorganEnvironment<'a> {
    topology: Cow<'a, TopologyBlock>,
    properties: Cow<'a, MoleculeProperties>,
}

impl PreparedMorganEnvironment<'_> {
    pub(crate) fn topology(&self) -> &TopologyBlock {
        &self.topology
    }

    pub(crate) fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }
}

pub(crate) fn prepare_morgan_environment<'a>(
    topology: &'a TopologyBlock,
    properties: &'a MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    include_chirality: bool,
) -> Result<PreparedMorganEnvironment<'a>, MorganError> {
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper chirality preparation
    // RDKit❗✔️:   const ROMol *lmol = &mol;
    // RDKit❗✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit❗✔️:   if (dp_fingerprintArguments->df_includeChirality &&
    // RDKit❗✔️:       !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit❗✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit❗✔️:     lmol = tmol.get();
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper chirality preparation
    // BEGIN RDKIT CPP FUNCTION ROMol::initFromOther relevant detached copy
    // RDKit❗✔️:   for (const auto oatom : other.atoms()) {
    // RDKit❗✔️:     constexpr bool updateLabel = false;
    // RDKit❗✔️:     constexpr bool takeOwnership = true;
    // RDKit❗✔️:     addAtom(oatom->copy(), updateLabel, takeOwnership);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (const auto obond : other.bonds()) {
    // RDKit❗✔️:     addBond(obond->copy(), true);
    // RDKit❗✔️:   }
    // RDKit❗✔️:     d_props = other.d_props;
    // END RDKIT CPP FUNCTION ROMol::initFromOther relevant detached copy
    // BEGIN RDKIT CPP FUNCTION MolOps::assignStereochemistry default property effect
    // RDKit❗✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // END RDKIT CPP FUNCTION MolOps::assignStereochemistry default property effect
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT void assignStereochemistry(
    // RDKit✔️✔️:     ROMol &mol, bool cleanIt = false, bool force = false,
    // RDKit✔️✔️:     bool flagPossibleStereoCenters = false);
    // Behavior review: test `_StereochemDone` by property presence only, so
    // both present value "0" and value "1" borrow the original record. The
    // missing-property branch assigns only to owned copies; the core call is
    // exactly cleanIt=false/flagPossibleStereoCenters=false and its typed
    // error is preserved. The source's computed done-property write happens
    // only after assignment succeeds. The fixed legacy profile is used; the
    // independent modern CIPLabeler branch is not substituted here.
    // Core's current false/false owner writes `_ChiralityPossible` for legal
    // centers on the private copy even though the pinned false branch does
    // not; Morgan's environment does not read this property. Preserve that
    // known owner delta explicitly instead of adding local chemistry cleanup.
    // Complexity review: the no-copy branches borrow both blocks. The source
    // deep-copies a full ROMol (including conformers and other collections);
    // this detached view clones the topology and molecule properties once
    // only when the source condition requires it. Coordinates are not read by
    // this assignment/environment path, so omitting their clone preserves
    // the Morgan result while reducing temporary allocation.
    if !include_chirality || properties.prop("_StereochemDone").is_some() {
        return Ok(PreparedMorganEnvironment {
            topology: Cow::Borrowed(topology),
            properties: Cow::Borrowed(properties),
        });
    }

    let assigned_topology =
        assign_legacy_stereochemistry_for_depiction(topology.clone(), valence, rings)?;
    let mut assigned_properties = properties.clone();
    assigned_properties.set_computed_prop("_StereochemDone", "1")?;

    Ok(PreparedMorganEnvironment {
        topology: Cow::Owned(assigned_topology),
        properties: Cow::Owned(assigned_properties),
    })
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::{PreparedMorganEnvironment, prepare_morgan_environment};
    use cosmolkit_core::{ValenceParams, assign_valence, fast_find_rings};
    use cosmolkit_model::{BondOrder, BondStereo, ChiralTag, TopologyBlock};
    use cosmolkit_smiles::{SmilesParseParams, SmilesRecord, parse_smiles};

    fn fixed_record(smiles: &str) -> SmilesRecord {
        parse_smiles(smiles, &SmilesParseParams::default()).expect("fixed I07 raw SMILES parses")
    }

    fn prepare<'a>(
        record: &'a SmilesRecord,
        include_chirality: bool,
    ) -> PreparedMorganEnvironment<'a> {
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed I07 input has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed I07 input has ring info");
        prepare_morgan_environment(
            &record.topology,
            &record.properties,
            &valence,
            &rings,
            include_chirality,
        )
        .expect("fixed I07 chirality preparation succeeds")
    }

    fn cip_code(topology: &TopologyBlock, atom_index: usize) -> Option<&str> {
        topology.atoms[atom_index]
            .prop("_CIPCode")
            .map(|value| value.as_string().expect("fixed CIP label is a string"))
    }

    #[test]
    fn fingerprint_morgan_i07_preparation_complete_source_product() {
        const DONE_VALUES: [Option<&str>; 3] = [None, Some("0"), Some("1")];
        const INPUT_CIP_LABELS: [Option<&str>; 4] = [None, Some("R"), Some("S"), Some("invalid")];

        let mut cases = 0;
        for include_chirality in [false, true] {
            for done_value in DONE_VALUES {
                for input_cip in INPUT_CIP_LABELS {
                    let mut record = fixed_record("F[C@H](Cl)Br");
                    if let Some(done_value) = done_value {
                        record
                            .properties
                            .set_prop("_StereochemDone", done_value)
                            .unwrap();
                    }
                    if let Some(input_cip) = input_cip {
                        record.topology.atoms[1]
                            .set_prop("_CIPCode", input_cip)
                            .unwrap();
                    }
                    let original_topology = record.topology.clone();
                    let original_properties = record.properties.clone();

                    let prepared = prepare(&record, include_chirality);
                    let source_assigns = include_chirality && done_value.is_none();
                    assert_eq!(
                        matches!(&prepared.topology, Cow::Owned(_)),
                        source_assigns,
                        "include_chirality={include_chirality}, done={done_value:?}, cip={input_cip:?}: topology ownership"
                    );
                    assert_eq!(
                        matches!(&prepared.properties, Cow::Owned(_)),
                        source_assigns,
                        "include_chirality={include_chirality}, done={done_value:?}, cip={input_cip:?}: property ownership"
                    );
                    let expected_cip = if source_assigns {
                        input_cip.or(Some("R"))
                    } else {
                        input_cip
                    };
                    assert_eq!(
                        prepared.topology.atoms[1]
                            .prop("_CIPCode")
                            .map(|value| value.as_string().unwrap()),
                        expected_cip,
                        "include_chirality={include_chirality}, done={done_value:?}, cip={input_cip:?}: CIP"
                    );
                    if source_assigns {
                        assert_eq!(prepared.properties().prop("_StereochemDone"), Some("1"));
                        assert!(prepared.properties().is_prop_computed("_StereochemDone"));
                    } else {
                        assert_eq!(prepared.properties().prop("_StereochemDone"), done_value);
                        assert!(!prepared.properties().is_prop_computed("_StereochemDone"));
                    }
                    assert_eq!(record.topology, original_topology);
                    assert_eq!(record.properties, original_properties);
                    cases += 1;
                }
            }
        }
        assert_eq!(cases, 24, "2 × 3 × 4 complete source-state product");

        for (smiles, expected_cip) in [("F[C@H](Cl)Br", "R"), ("F[C@@H](Cl)Br", "S")] {
            let record = fixed_record(smiles);
            let original = record.topology.clone();
            let prepared = prepare(&record, true);
            assert!(matches!(&prepared.topology, Cow::Owned(_)));
            assert_eq!(
                prepared.topology.atoms[1]
                    .prop("_CIPCode")
                    .map(|value| value.as_string().unwrap()),
                Some(expected_cip),
                "pinned legacy label for {smiles}"
            );
            assert_eq!(record.topology, original);
        }

        let invalid = fixed_record("F[C@](Cl)(Cl)Br");
        let invalid_tag = invalid.topology.atoms[1].chiral_tag();
        let invalid_prepared = prepare(&invalid, true);
        assert_eq!(invalid_tag, ChiralTag::TetrahedralCcw);
        assert_eq!(
            invalid_prepared.topology.atoms[1].chiral_tag(),
            ChiralTag::TetrahedralCcw,
            "cleanIt=false retains the specified tag on a duplicate-ligand center"
        );
        assert_eq!(cip_code(invalid_prepared.topology.as_ref(), 1), None);

        let unspecified = fixed_record("FC(Cl)Br");
        let unspecified_prepared = prepare(&unspecified, true);
        assert_eq!(
            unspecified_prepared.topology.atoms[1].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(
            unspecified_prepared.topology.atoms[1].prop("_CIPCode"),
            None
        );

        for (smiles, expected) in [("F/C=C/F", BondStereo::E), (r"F/C=C\F", BondStereo::Z)] {
            let record = fixed_record(smiles);
            let original = record.topology.clone();
            let prepared = prepare(&record, true);
            let double_bond = prepared
                .topology
                .bonds
                .iter()
                .find(|bond| bond.order() == BondOrder::Double)
                .expect("fixed E/Z input has a double bond");
            assert_eq!(
                double_bond.stereo(),
                expected,
                "source direction assignment: {smiles}"
            );
            assert_eq!(record.topology, original);
        }
    }
}
