//! InChI facade: checked model transport and calls to the single InChI owner.
mod toolkit;

use crate::{InchiError, InchiErrorKind, Molecule};

/// Controls the source-defined InChI post-processing stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InchiReadParams {
    pub sanitize: bool,
    pub remove_hs: bool,
}

impl Default for InchiReadParams {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hs: true,
        }
    }
}

impl InchiReadParams {
    pub const fn new(sanitize: bool, remove_hs: bool) -> Self {
        Self {
            sanitize,
            remove_hs,
        }
    }
    pub const fn sanitize(&self) -> bool {
        self.sanitize
    }
    pub const fn remove_hs(&self) -> bool {
        self.remove_hs
    }
    pub fn set_sanitize(&mut self, value: bool) {
        self.sanitize = value;
    }
    pub fn set_remove_hs(&mut self, value: bool) {
        self.remove_hs = value;
    }
}

/// Official InChI option string; an empty string selects the engine defaults.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct InchiWriteParams {
    pub options: String,
}

impl InchiWriteParams {
    pub fn new(options: String) -> Self {
        Self { options }
    }
    pub fn options(&self) -> &str {
        &self.options
    }
    pub fn set_options(&mut self, value: String) {
        self.options = value;
    }
}

fn text(bytes: Vec<u8>, operation: &'static str) -> Result<String, InchiError> {
    if bytes.is_empty() {
        return Err(InchiError {
            operation,
            kind: InchiErrorKind::InvalidInput,
            detail: "the InChI engine returned no identifier".into(),
        });
    }
    String::from_utf8(bytes).map_err(|error| InchiError {
        operation,
        kind: InchiErrorKind::InvalidSourceOutput,
        detail: error.to_string(),
    })
}

/// Converts an InChI string to its InChIKey without constructing a molecule.
pub fn inchi_to_key(inchi: &str) -> Result<String, InchiError> {
    text(
        cosmolkit_inchi::inchi_to_inchi_key(inchi.as_bytes())?.key,
        "inchi_to_key",
    )
}

impl Molecule {
    /// Parses InChI using sanitization and hydrogen removal by default.
    pub fn from_inchi(text: &str) -> Result<Self, InchiError> {
        Self::from_inchi_with_params(text, &InchiReadParams::default())
    }

    pub fn from_inchi_with_params(
        text: &str,
        params: &InchiReadParams,
    ) -> Result<Self, InchiError> {
        let mut toolkit = toolkit::Toolkit::default();
        let output = cosmolkit_inchi::mol_from_inchi(
            &mut toolkit,
            text.as_bytes(),
            params.sanitize,
            params.remove_hs,
        )?;
        let graph = output.molecule.ok_or_else(|| InchiError {
            operation: "from_inchi",
            kind: InchiErrorKind::InvalidInput,
            detail: String::from_utf8_lossy(&output.return_values.message).into_owned(),
        })?;
        let (_, coordinates) = toolkit::model(&graph).map_err(toolkit::public_error)?;
        let topology = toolkit.final_topology.ok_or_else(|| InchiError {
            operation: "from_inchi",
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: "InChI parsing did not execute final stereo assignment".into(),
        })?;
        Self::from_parsed_parts_with_derived_state(
            topology,
            coordinates,
            Default::default(),
            toolkit::cached_valence(&graph),
            toolkit.rings,
        )
        .map_err(|error| InchiError {
            operation: "from_inchi",
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: error.to_string(),
        })
    }

    /// Generates InChI without changing the receiver or its caches.
    pub fn to_inchi(&self) -> Result<String, InchiError> {
        self.to_inchi_with_params(&InchiWriteParams::default())
    }

    pub fn to_inchi_with_params(&self, params: &InchiWriteParams) -> Result<String, InchiError> {
        let graph = toolkit::graph(
            self.topology(),
            self.coordinate_block_runtime(),
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(toolkit::public_error)?;
        let output = cosmolkit_inchi::mol_to_inchi(
            &mut toolkit::Toolkit::default(),
            &graph,
            Some(params.options.as_bytes()),
        )?;
        text(output.inchi, "to_inchi")
    }

    /// Generates an InChIKey without changing the receiver.
    pub fn to_inchi_key(&self) -> Result<String, InchiError> {
        self.to_inchi_key_with_params(&InchiWriteParams::default())
    }

    pub fn to_inchi_key_with_params(
        &self,
        params: &InchiWriteParams,
    ) -> Result<String, InchiError> {
        let graph = toolkit::graph(
            self.topology(),
            self.coordinate_block_runtime(),
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(toolkit::public_error)?;
        let output = cosmolkit_inchi::mol_to_inchi_key(
            &mut toolkit::Toolkit::default(),
            &graph,
            Some(params.options.as_bytes()),
        )?;
        text(output.key, "to_inchi_key")
    }
}

#[cfg(test)]
mod tests {
    use crate::{AtomSpec, BondOrder, BondSpec, BondStereo, Element, MoleculeBuilder};
    #[test]
    fn rdkit_2026_03_inchi_unknown_export_core_boundary() {
        // Explicit MolToInchi precondition, independent of CXSMILES parsing.
        // .6 inchi.cpp:2037–2078 consumes Any and adjacency, not stereoAtoms.
        // The separate nine official SMILES regression remains unchanged.
        for reverse in [false, true] {
            let mut builder = MoleculeBuilder::new();
            let fluorine = builder.add_atom(AtomSpec::new(Element::F));
            let left = builder.add_atom(AtomSpec::new(Element::C));
            let right = builder.add_atom(AtomSpec::new(Element::C));
            let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
            builder
                .add_bond(BondSpec::new(fluorine, left, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(
                    BondSpec::new(
                        if reverse { right } else { left },
                        if reverse { left } else { right },
                        BondOrder::Double,
                    )
                    .with_stereo(BondStereo::Any),
                )
                .unwrap();
            builder
                .add_bond(BondSpec::new(right, chlorine, BondOrder::Single))
                .unwrap();
            let molecule = builder.build().unwrap();
            assert_eq!(molecule.bonds()[1].stereo(), BondStereo::Any);
            assert_eq!(molecule.bonds()[1].stereo_atoms(), None);
            let before = molecule.clone();
            let generated = molecule
                .to_inchi_with_params(&InchiWriteParams::new("-SUU".into()))
                .unwrap();
            assert!(!generated.is_empty());
            assert_eq!(molecule, before);
            let restored = Molecule::from_inchi(&generated).unwrap();
            let unknown = restored
                .bonds()
                .iter()
                .filter(|bond| bond.stereo() == BondStereo::Any)
                .collect::<Vec<_>>();
            assert_eq!(unknown.len(), 1);
            assert_eq!(unknown[0].order(), BondOrder::Double);
            let controllers = unknown[0]
                .stereo_atoms()
                .expect("INCHI01 UNKNOWN import retains both controllers");
            assert_ne!(controllers[0], controllers[1]);
            assert!(
                controllers
                    .iter()
                    .all(|atom| atom.index() < restored.num_atoms())
            );
        }
    }

    #[cfg(feature = "cap-smiles")]
    fn check_rdkit_2026_03_inchi_unknown_roundtrip(smiles: &str) {
        let molecule = Molecule::from_smiles(smiles).unwrap();
        assert!(
            molecule
                .bonds()
                .iter()
                .any(|b| b.stereo() == BondStereo::Any),
            "{smiles}"
        );
        // RDKit .6 inchi.cpp:1751 copies the input molecule before modifying it.
        // Compare the complete modeled value: the current detached SMILES writer
        // independently rejects Any stereo, so serialization cannot observe this input.
        let before = molecule.clone();
        let generated = molecule
            .to_inchi_with_params(&InchiWriteParams::new("-SUU".into()))
            .unwrap();
        assert!(!generated.is_empty());
        assert_eq!(molecule, before);
        let restored = Molecule::from_inchi(&generated).unwrap();
        let unknown = restored
            .bonds()
            .iter()
            .filter(|b| b.stereo() == BondStereo::Any)
            .collect::<Vec<_>>();
        assert!(!unknown.is_empty(), "{smiles}: {generated:?}");
        assert!(
            unknown.iter().all(|b| b.stereo_atoms().is_some()),
            "{smiles}"
        );
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_01() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("CSC1=NSC(CC=NC2=CC=CC=C2)=C1C#N |w:8.8|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_02() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("O/N=C/c1ccccc1 |w:1.1|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_03() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("OC(=O)/C=C/c1ccccc1 |w:3.3|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_04() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("O=C(/C=C/c1ccccc1)c1ccccc1 |w:2.2|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_05() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("C/C=C/C=O |w:1.1|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_06() {
        check_rdkit_2026_03_inchi_unknown_roundtrip(
            "CC/C(=C(/c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1 |w:2.2|",
        );
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_07() {
        check_rdkit_2026_03_inchi_unknown_roundtrip(
            "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=O)C(C)(C)CCC1 |w:3.3|",
        );
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_08() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("O=C(O)[C@@H](CC=Cc1ccccc1)N |w:4.4|");
    }
    #[cfg(feature = "cap-smiles")]
    #[test]
    fn rdkit_2026_03_inchi_unknown_roundtrip_upstream_09() {
        check_rdkit_2026_03_inchi_unknown_roundtrip("C[C@H](O)/C(=N/O)c1ccccc1 |w:3.3|");
    }
    use super::*;

    #[test]
    fn inchi_feature_alone_parses_and_generates_identifiers() {
        let identifier = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3";
        let molecule = Molecule::from_inchi(identifier).unwrap();
        assert_eq!(molecule.to_inchi().unwrap(), identifier);
        assert_eq!(
            molecule.to_inchi_key().unwrap(),
            "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        );
        assert!(Molecule::from_inchi("not InChI").is_err());
        assert!(inchi_to_key("not InChI").is_err());
    }

    #[cfg(feature = "cap-smiles")]
    #[test]
    fn inchi_identifiers_round_trip_and_preserve_the_receiver() {
        // Fixed identifiers cross-checked with RDKit 2026.03.1 InChI.
        for (smiles, expected, key) in [
            (
                "CCO",
                "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            ),
            (
                "c1ccccc1",
                "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
                "UHOVQNZJYSORNB-UHFFFAOYSA-N",
            ),
            (
                "[13CH4]",
                "InChI=1S/CH4/h1H4/i1+1",
                "VNWKTOKETHGBQD-OUBTZVSYSA-N",
            ),
            (
                "F[C@H](Cl)Br",
                "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1",
                "YNKZSBSRKWVMEZ-SFOWXEAESA-N",
            ),
            (
                "F/C=C/F",
                "InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+",
                "WFLOTYSKFUPZQB-OWOJBTEDSA-N",
            ),
            (
                "[Na+].[Cl-]",
                "InChI=1S/ClH.Na/h1H;/q;+1/p-1",
                "FAPWRFPIFSIZLT-UHFFFAOYSA-M",
            ),
        ] {
            let molecule = Molecule::from_smiles(smiles).unwrap();
            let before = molecule.to_smiles().unwrap();
            assert_eq!(molecule.to_inchi().unwrap(), expected, "{smiles}");
            assert_eq!(molecule.to_inchi_key().unwrap(), key, "{smiles}");
            assert_eq!(inchi_to_key(expected).unwrap(), key);
            let parsed = Molecule::from_inchi(expected).unwrap();
            assert_eq!(parsed.to_inchi().unwrap(), expected, "{smiles}");
            assert_eq!(molecule.to_smiles().unwrap(), before);
        }
    }
}
