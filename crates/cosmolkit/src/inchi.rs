//! InChI facade: checked model transport and calls to the single InChI owner.
mod toolkit;

use crate::{InchiError, InchiErrorKind, Molecule};

/// Controls the source-defined InChI post-processing stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InchiReadParams {
    pub sanitize: bool,
    pub remove_hydrogens: bool,
}

impl Default for InchiReadParams {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hydrogens: true,
        }
    }
}

impl InchiReadParams {
    pub const fn new(sanitize: bool, remove_hydrogens: bool) -> Self {
        Self {
            sanitize,
            remove_hydrogens,
        }
    }
    pub const fn sanitize(&self) -> bool {
        self.sanitize
    }
    pub const fn remove_hydrogens(&self) -> bool {
        self.remove_hydrogens
    }
    pub fn set_sanitize(&mut self, value: bool) {
        self.sanitize = value;
    }
    pub fn set_remove_hydrogens(&mut self, value: bool) {
        self.remove_hydrogens = value;
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
pub fn inchi_to_inchi_key(inchi: &str) -> Result<String, InchiError> {
    text(
        cosmolkit_inchi::inchi_to_inchi_key(inchi.as_bytes())?.key,
        "inchi_to_inchi_key",
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
            params.remove_hydrogens,
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
        assert!(inchi_to_inchi_key("not InChI").is_err());
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
            assert_eq!(inchi_to_inchi_key(expected).unwrap(), key);
            let parsed = Molecule::from_inchi(expected).unwrap();
            assert_eq!(parsed.to_inchi().unwrap(), expected, "{smiles}");
            assert_eq!(molecule.to_smiles().unwrap(), before);
        }
    }
}
