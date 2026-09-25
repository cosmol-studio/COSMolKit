//! Pinned RDKit ring-system coordinate templates over canonical query graphs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cosmolkit_core::{
    PathError, RingFindingError, RingInfo, connected_components,
    symmetrize_sssr_with_options_from_parts,
};
use cosmolkit_model::{QueryGraph, TopologyBlock, TopologyValidationError};
use cosmolkit_search::{SmartsParseError, SmartsParseParams, parse_smarts};

const DEFAULT_TEMPLATE_ROWS: &str = include_str!("../assets/default_templates.cxsmarts");

#[derive(Debug)]
pub(crate) enum TemplateError {
    InvalidDefaultRow {
        index: usize,
        source: SmartsParseError,
    },
    InvalidTopology {
        index: usize,
        source: TopologyValidationError,
    },
    RingInitialization {
        index: usize,
        source: RingFindingError,
    },
    ConnectedComponents {
        index: usize,
        source: PathError,
    },
    ExternalOpen {
        path: String,
        source: std::io::Error,
    },
    ExternalRead {
        path: String,
        source: std::io::Error,
    },
    ExternalInvalidSmarts {
        path: String,
        line: usize,
        source: SmartsParseError,
    },
    MissingCoordinates {
        row: String,
    },
    ThreeDimensionalCoordinates {
        row: String,
    },
    MultipleFragments {
        row: String,
    },
    NotRingSystem {
        row: String,
    },
}

#[derive(Debug)]
pub(crate) struct RingTemplate {
    pub(crate) query: QueryGraph,
    pub(crate) rings: RingInfo,
}

#[derive(Debug, Default)]
pub(crate) struct CoordinateTemplates {
    buckets: BTreeMap<usize, Vec<RingTemplate>>,
}

impl CoordinateTemplates {
    pub(crate) fn new_default() -> Result<Self, TemplateError> {
        let mut templates = Self::default();
        templates.load_default_templates()?;
        Ok(templates)
    }

    fn clear_templates(&mut self) {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::clearTemplates
        // RDKit❗❌: void clearTemplates() {
        // RDKit❗❌:   for (auto &[atom_cout, romols] : m_templates) {
        // RDKit❗❌:     romols.clear();
        // RDKit❗❌:   }
        // RDKit❗❌:   m_templates.clear();
        // RDKit❗❌: }
        // Behavior: clearing the owning map drops all bucket contents; no
        // external references to a mutable singleton are exposed here.
        // Complexity: BTreeMap teardown and individual QueryGraph drops have
        // higher pointer-chasing cost than the source unordered map.
        self.buckets.clear();
    }

    pub(crate) fn load_default_templates(&mut self) -> Result<(), TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::loadDefaultTemplates
        // RDKit❗❌: void loadDefaultTemplates() {
        // RDKit❗❌:   clearTemplates();
        self.clear_templates();
        self.load_default_rows(DEFAULT_TEMPLATE_ROWS)
    }

    // The row parameter lets the owner test the same failure path without
    // modifying the pinned, generated production corpus.
    fn load_default_rows(&mut self, rows: &str) -> Result<(), TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::loadDefaultTemplates (row loop)
        // RDKit❗❌:   // load default templates into m_templates map by atom count
        // RDKit❗❌:   for (const auto &smarts : TEMPLATE_SMARTS) {
        // RDKit❗❌:     std::shared_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
        // RDKit❗❌:     if (!mol) {
        // RDKit❗❌:       continue;
        // RDKit❗❌:     }
        // RDKit❗❌:     // Initialize ring info using symmetrizeSSSR to match depictor ring counting
        // RDKit❗❌:     RDKit::VECT_INT_VECT arings;
        // RDKit❗❌:     RDKit::MolOps::symmetrizeSSSR(*mol, arings);
        // RDKit❗❌:     m_templates[mol->getNumAtoms()].push_back(mol);
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // Behavior: source-null SMARTS rows are skipped without disturbing
        // earlier or later rows. All 578 pinned default rows parse in both
        // fixed RDKit 2026.03.1 and this parser. A new parse failure in that
        // known source-valid corpus is a parser regression, not source null.
        // Arbitrary SMARTS parser-error equivalence is not established here.
        // Complexity: a temporary validated carrier topology and BTreeMap
        // buckets add O(V+E) copying per row beyond RDKit's in-place RWMol.
        for (index, text) in rows.lines().enumerate() {
            let (template, _) = match Self::parse_template_row(text, index) {
                Ok(row) => row,
                Err(TemplateError::InvalidDefaultRow { .. }) if rows != DEFAULT_TEMPLATE_ROWS => {
                    continue;
                }
                Err(error) => return Err(error),
            };
            self.buckets
                .entry(template.query.num_atoms())
                .or_default()
                .push(template);
        }
        Ok(())
    }

    pub(crate) fn templates_of_size(&self, atom_count: usize) -> Option<&[RingTemplate]> {
        self.buckets.get(&atom_count).map(Vec::as_slice)
    }

    pub(crate) fn has_template_of_size(&self, atom_count: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::hasTemplateOfSize
        // RDKit✔️❌: bool hasTemplateOfSize(unsigned int atomCount) {
        // RDKit✔️❌:   if (m_templates.find(atomCount) != m_templates.end()) {
        // RDKit✔️❌:     return true;
        // RDKit✔️❌:   }
        // RDKit✔️❌:   return false;
        // RDKit✔️❌: }
        // Behavior: bucket existence, including a bucket made by a failed-size
        // lookup, is distinct from whether the bucket contains templates.
        // Complexity: BTreeMap lookup is O(log n) versus unordered-map O(1)
        // average; this private registry has a bounded number of atom sizes.
        self.buckets.contains_key(&atom_count)
    }

    pub(crate) fn matching_templates(&mut self, atom_count: usize) -> &[RingTemplate] {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::getMatchingTemplates
        // RDKit✔️❌: const std::vector<std::shared_ptr<RDKit::ROMol>> &getMatchingTemplates(
        // RDKit✔️❌:     unsigned int atomCount) {
        // RDKit✔️❌:   return m_templates[atomCount];
        // RDKit✔️❌: }
        // Behavior: a missing lookup installs an empty bucket, affecting
        // subsequent hasTemplateOfSize just as source operator[] does.
        // Complexity: BTreeMap entry is O(log n), unlike average O(1) source.
        self.buckets.entry(atom_count).or_default().as_slice()
    }

    pub(crate) fn set_ring_system_templates(&mut self, path: &Path) -> Result<(), TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::setRingSystemTemplates
        // RDKit✔️❌: void CoordinateTemplates::setRingSystemTemplates(
        // RDKit✔️❌:     const std::string &templatePath) {
        // RDKit✔️❌:   // Try loading templates in from this directory, if unsuccessful, keep current
        // RDKit✔️❌:   // templates
        // RDKit✔️❌:   std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
        // RDKit✔️❌:       templates;
        // RDKit✔️❌:   loadTemplatesFromPath(templatePath, templates);
        // RDKit✔️❌:   clearTemplates();
        // RDKit✔️❌:   m_templates = std::move(templates);
        // RDKit✔️❌: }
        // Behavior: no mutation occurs until the entire external file loads.
        // Complexity: BTreeMap buckets have O(log n) insertion rather than
        // average O(1); map replacement is a move in both implementations.
        let staged = Self::load_templates_from_path(path)?;
        self.clear_templates();
        self.buckets = staged;
        Ok(())
    }

    pub(crate) fn add_ring_system_templates(&mut self, path: &Path) -> Result<(), TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::addRingSystemTemplates
        // RDKit✔️❌: void CoordinateTemplates::addRingSystemTemplates(
        // RDKit✔️❌:     const std::string &templatePath) {
        // RDKit✔️❌:   // Try loading templates in from this directory, if unsuccessful, keep current
        // RDKit✔️❌:   // templates
        // RDKit✔️❌:   std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
        // RDKit✔️❌:       templates;
        // RDKit✔️❌:   loadTemplatesFromPath(templatePath, templates);
        // RDKit✔️❌:   for (auto &kv : templates) {
        // RDKit✔️❌:     m_templates[kv.first].insert(m_templates[kv.first].begin(),
        // RDKit✔️❌:                                  kv.second.begin(), kv.second.end());
        // RDKit✔️❌:   }
        // RDKit✔️❌: }
        // Behavior: each staged size group prepends in its original file order;
        // an invalid file leaves every existing bucket untouched.
        // Complexity: per-bucket front insertion shifts existing entries, as
        // in source vector::insert; BTreeMap lookup is logarithmic.
        let staged = Self::load_templates_from_path(path)?;
        for (size, rows) in staged {
            self.buckets.entry(size).or_default().splice(0..0, rows);
        }
        Ok(())
    }

    pub(crate) fn template_count(&self) -> usize {
        self.buckets.values().map(Vec::len).sum()
    }

    fn parse_template_row(
        text: &str,
        index: usize,
    ) -> Result<(RingTemplate, TopologyBlock), TemplateError> {
        // RDKit❗❌:     std::shared_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
        // RDKit❗❌:     // Initialize ring info using symmetrizeSSSR to match depictor ring counting
        // RDKit❗❌:     RDKit::VECT_INT_VECT arings;
        // RDKit❗❌:     RDKit::MolOps::symmetrizeSSSR(*mol, arings);
        // Behavior: parse failure is returned to the caller for the source
        // null-row rule. The fixed corpus has no null parses in either pinned
        // reference environment or this parser; equivalence for arbitrary
        // external SMARTS text is not established by that corpus check.
        // Complexity: the canonical QueryGraph is retained, while a validated
        // carrier topology is copied once for core structural algorithms.
        let query = parse_smarts(text, &SmartsParseParams::default())
            .map_err(|source| TemplateError::InvalidDefaultRow { index, source })?;
        let atoms = query
            .atoms()
            .iter()
            .map(|atom| atom.atom().clone())
            .collect();
        let bonds = query
            .bonds()
            .iter()
            .map(|bond| bond.bond().clone())
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .map_err(|source| TemplateError::InvalidTopology { index, source })?;
        let rings = symmetrize_sssr_with_options_from_parts(
            query.num_atoms(),
            &topology.bonds,
            &topology.adjacency,
            false,
            false,
        )
        .map_err(|source| TemplateError::RingInitialization { index, source })?;
        Ok((RingTemplate { query, rings }, topology))
    }

    fn assert_valid_template(
        template: &RingTemplate,
        topology: &TopologyBlock,
        row: &str,
    ) -> Result<(), TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::assertValidTemplate
        // RDKit❗❌: void CoordinateTemplates::assertValidTemplate(RDKit::ROMol &mol,
        // RDKit❗❌:                                               const std::string &smiles) {
        // RDKit❗❌:   // template must have 2D coordinates
        // RDKit❗❌:   if (mol.getNumConformers() == 0) {
        // RDKit❗❌:     std::string msg = "Template missing coordinates: " + smiles;
        // RDKit❗❌:     throw RDDepict::DepictException(msg);
        // RDKit❗❌:   }
        let query = &template.query;
        if query.coordinates_2d().is_none() && query.conformers_3d().is_empty() {
            return Err(TemplateError::MissingCoordinates { row: row.into() });
        }
        // RDKit❗❌:   if (mol.getConformer().is3D()) {
        // RDKit❗❌:     std::string msg =
        // RDKit❗❌:         "Template has 3D coordinates, 2D coordinates required: " + smiles;
        // RDKit❗❌:     throw RDDepict::DepictException(msg);
        // RDKit❗❌:   }
        if query.coordinates_2d().is_none()
            && query
                .conformers_3d()
                .first()
                .is_some_and(|conformer| conformer.is_3d())
        {
            return Err(TemplateError::ThreeDimensionalCoordinates { row: row.into() });
        }
        // RDKit❗❌:   // Make sure this template is a single ring system (spiro'd ring systems are
        // RDKit❗❌:   // OK). We can check this by ensuring that every bond is in a ring and that
        // RDKit❗❌:   // there is only one connected component in the molecular graph.
        // RDKit❗❌:   if (RDKit::MolOps::getMolFrags(mol).size() != 1) {
        // RDKit❗❌:     std::string msg =
        // RDKit❗❌:         "Template consists of multiple fragments, single fragment required: " +
        // RDKit❗❌:         smiles;
        // RDKit❗❌:     throw RDDepict::DepictException(msg);
        // RDKit❗❌:   }
        let components = connected_components(topology)
            .map_err(|source| TemplateError::ConnectedComponents { index: 0, source })?;
        if components.components.len() != 1 {
            return Err(TemplateError::MultipleFragments { row: row.into() });
        }
        // RDKit❗❌:   if (mol.getNumAtoms() == 1) {
        // RDKit❗❌:     std::string msg = "Template is not a ring system: " + smiles;
        // RDKit❗❌:     throw RDDepict::DepictException(msg);
        // RDKit❗❌:   }
        if query.num_atoms() == 1 {
            return Err(TemplateError::NotRingSystem { row: row.into() });
        }
        // RDKit❗❌:   auto ri = mol.getRingInfo();
        // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
        // RDKit❗❌:     if (!ri->numBondRings(i)) {
        // RDKit❗❌:       std::string msg = "Template is not a ring system: " + smiles;
        // RDKit❗❌:       throw RDDepict::DepictException(msg);
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // Behavior: uses source check order over the authoritative query's
        // carrier structure; no predicate or coordinate projection occurs.
        // Complexity: core components/ring lookup has comparable graph scale,
        // but temporary topology copying remains slower than source ROMol.
        if topology
            .bonds
            .iter()
            .any(|bond| template.rings.num_bond_rings(bond.id()) == 0)
        {
            return Err(TemplateError::NotRingSystem { row: row.into() });
        }
        Ok(())
    }

    pub(crate) fn load_templates_from_path(
        path: &Path,
    ) -> Result<BTreeMap<usize, Vec<RingTemplate>>, TemplateError> {
        // BEGIN RDKIT CPP FUNCTION CoordinateTemplates::loadTemplatesFromPath
        // RDKit❗❌: void CoordinateTemplates::loadTemplatesFromPath(
        // RDKit❗❌:     const std::string &templatePath,
        // RDKit❗❌:     std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
        // RDKit❗❌:         &templates) {
        // RDKit❗❌:   std::ifstream cxsmiles(templatePath);
        // RDKit❗❌:   if (!cxsmiles) {
        // RDKit❗❌:     std::string msg = "Could not open file " + templatePath;
        // RDKit❗❌:     throw RDDepict::DepictException(msg);
        // RDKit❗❌:   }
        let file = File::open(path).map_err(|source| TemplateError::ExternalOpen {
            path: path.display().to_string(),
            source,
        })?;
        let mut reader = BufReader::new(file);
        let mut templates = BTreeMap::<usize, Vec<RingTemplate>>::new();
        let mut row = String::new();
        let mut line = 0;
        // RDKit❗❌:   // Try loading templates in from this directory, if unsuccessful, keep current
        // RDKit❗❌:   // templates
        // RDKit❗❌:   std::string line;
        // RDKit❗❌:   while (std::getline(cxsmiles, line)) {
        loop {
            row.clear();
            let read =
                reader
                    .read_line(&mut row)
                    .map_err(|source| TemplateError::ExternalRead {
                        path: path.display().to_string(),
                        source,
                    })?;
            if read == 0 {
                break;
            }
            if row.ends_with('\n') {
                row.pop();
            }
            line += 1;
            // RDKit❗❌:     RDKit::ROMol *mol_ptr = RDKit::SmartsToMol(line);
            // RDKit❗❌:     if (!mol_ptr) {
            // RDKit❗❌:       std::string msg =
            // RDKit❗❌:           "Could not load templates from " + templatePath + ": Invalid smarts";
            // RDKit❗❌:       cxsmiles.close();
            // RDKit❗❌:       throw RDDepict::DepictException(msg);
            // RDKit❗❌:     }
            let (template, topology) =
                Self::parse_template_row(&row, line - 1).map_err(|error| match error {
                    TemplateError::InvalidDefaultRow { source, .. } => {
                        TemplateError::ExternalInvalidSmarts {
                            path: path.display().to_string(),
                            line,
                            source,
                        }
                    }
                    other => other,
                })?;
            // RDKit❗❌:     std::shared_ptr<RDKit::ROMol> mol(mol_ptr);
            // RDKit❗❌:     // Initialize ring info using symmetrizeSSSR to match depictor ring counting
            // RDKit❗❌:     RDKit::VECT_INT_VECT arings;
            // RDKit❗❌:     RDKit::MolOps::symmetrizeSSSR(*mol, arings);
            // RDKit❗❌:     try {
            // RDKit❗❌:       assertValidTemplate(*mol, line);
            // RDKit❗❌:     } catch (RDDepict::DepictException &e) {
            // RDKit❗❌:       cxsmiles.close();
            // RDKit❗❌:       throw e;
            // RDKit❗❌:     }
            Self::assert_valid_template(&template, &topology, &row)?;
            // RDKit❗❌:     templates[mol->getNumAtoms()].push_back(mol);
            templates
                .entry(template.query.num_atoms())
                .or_default()
                .push(template);
        }
        // RDKit❗❌:   }
        // RDKit❗❌:   cxsmiles.close();
        // RDKit❗❌: }
        // Behavior: stage in a private map; the caller's registry is never
        // touched if opening, reading, parsing or validation fails.
        // Complexity: line buffering is comparable; temporary carrier copy
        // and BTreeMap insertion cost more than source ROMol/unordered_map.
        Ok(templates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT_TEST_PATH: AtomicUsize = AtomicUsize::new(0);

    fn with_external_rows<T>(rows: &str, test: impl FnOnce(&Path) -> T) -> T {
        let index = NEXT_TEST_PATH.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-depict-templates-{}-{index}.cxsmarts",
            std::process::id()
        ));
        std::fs::write(&path, rows).expect("create owned temporary fixture");
        let result = test(&path);
        std::fs::remove_file(&path).expect("remove owned temporary fixture");
        result
    }

    #[test]
    fn default_templates_full_pinned_corpus_preserves_order_queries_and_coordinates() {
        let templates = CoordinateTemplates::new_default().expect("all pinned source rows parse");
        let rows: Vec<_> = DEFAULT_TEMPLATE_ROWS.lines().collect();
        assert_eq!(rows.len(), 578, "RDKit 2026.03.1 TemplateSmarts.h corpus");
        assert_eq!(templates.template_count(), rows.len());

        let mut per_size_index = BTreeMap::<usize, usize>::new();
        for (source_index, row) in rows.iter().enumerate() {
            let expected = parse_smarts(row, &SmartsParseParams::default())
                .unwrap_or_else(|error| panic!("source row {source_index}: {error}"));
            let next = per_size_index.entry(expected.num_atoms()).or_default();
            let stored = &templates.templates_of_size(expected.num_atoms()).unwrap()[*next];
            assert_eq!(stored.query, expected, "source row {source_index}");
            assert_eq!(stored.query.atoms(), expected.atoms());
            assert_eq!(stored.query.bonds(), expected.bonds());
            assert!(stored.rings.is_symm_sssr(), "source row {source_index}");
            let stored_coordinates = stored.query.coordinate_block(None);
            let expected_coordinates = expected.coordinate_block(None);
            assert_eq!(
                stored_coordinates.conformers_2d.len(),
                expected_coordinates.conformers_2d.len()
            );
            assert_eq!(
                stored_coordinates.conformers_3d.len(),
                expected_coordinates.conformers_3d.len()
            );
            for (actual, reference) in stored_coordinates
                .conformers_3d
                .iter()
                .zip(&expected_coordinates.conformers_3d)
            {
                assert_eq!(actual.is_3d(), reference.is_3d());
                for (actual_row, reference_row) in
                    actual.coordinates().iter().zip(reference.coordinates())
                {
                    assert_eq!(
                        actual_row.map(f64::to_bits),
                        reference_row.map(f64::to_bits)
                    );
                }
            }
            *next += 1;
        }
        assert_eq!(per_size_index.len(), templates.buckets.len());
        for (size, count) in per_size_index {
            assert_eq!(templates.templates_of_size(size).unwrap().len(), count);
        }
    }

    #[test]
    fn default_templates_source_null_parse_skips_only_that_row() {
        let mut templates = CoordinateTemplates::default();
        templates.load_default_rows("C1CC1\n[\nN1CC1\n").unwrap();
        assert_eq!(templates.template_count(), 2);
        assert_eq!(templates.templates_of_size(3).unwrap().len(), 2);
    }

    #[test]
    fn external_templates_valid_rows_keep_query_coordinates_and_spiro_ring() {
        // RDKit 2026.03.1 MolFromSmarts produces a non-3D conformer for each
        // row and 3 / 7 atoms respectively; Templates.cpp permits spiro rings.
        let triangle = "[C;R]1CC1 |(0,0,;1,0,;0.5,0.866,)|";
        let spiro = "C1CCC12CCC2 |(0,0,;1,0,;1,1,;0,1,;1,2,;2,2,;2,1,)|";
        with_external_rows(&format!("{triangle}\n{spiro}\n"), |path| {
            let loaded = CoordinateTemplates::load_templates_from_path(path).unwrap();
            assert_eq!(loaded.get(&3).unwrap().len(), 1);
            assert_eq!(loaded.get(&7).unwrap().len(), 1);
            let expected = parse_smarts(triangle, &SmartsParseParams::default()).unwrap();
            let actual = &loaded[&3][0].query;
            assert_eq!(actual, &expected);
            assert_eq!(
                actual.atom(0).unwrap().predicate(),
                expected.atom(0).unwrap().predicate()
            );
            assert!(loaded[&7][0].rings.num_rings() >= 2);
            assert!(!actual.conformers_3d()[0].is_3d());
            assert_eq!(
                actual.conformers_3d()[0].coordinates()[0].map(f64::to_bits),
                expected.conformers_3d()[0].coordinates()[0].map(f64::to_bits)
            );
        });
    }

    #[test]
    fn external_templates_source_validation_order_and_errors() {
        let cases = [
            ("[\n", "smarts"),
            ("C1CC1\n", "missing"),
            ("C1CC1 |(0,0,0;1,0,0;0.5,0.866,1)|\n", "three_d"),
            (
                "C1CC1.C1CC1 |(0,0,;1,0,;0.5,0.866,;3,0,;4,0,;3.5,0.866,)|\n",
                "fragments",
            ),
            ("C |(0,0,)|\n", "nonring"),
            ("C1CC1C |(0,0,;1,0,;0.5,0.866,;1.5,1.5,)|\n", "nonring"),
        ];
        for (row, expected) in cases {
            with_external_rows(row, |path| {
                let error = CoordinateTemplates::load_templates_from_path(path).unwrap_err();
                assert!(
                    match expected {
                        "smarts" =>
                            matches!(error, TemplateError::ExternalInvalidSmarts { line: 1, .. }),
                        "missing" => matches!(error, TemplateError::MissingCoordinates { .. }),
                        "three_d" =>
                            matches!(error, TemplateError::ThreeDimensionalCoordinates { .. }),
                        "fragments" => matches!(error, TemplateError::MultipleFragments { .. }),
                        "nonring" => matches!(error, TemplateError::NotRingSystem { .. }),
                        _ => false,
                    },
                    "case {expected}: {error:?}"
                );
            });
        }
        let missing = std::env::temp_dir().join(format!(
            "cosmolkit-depict-missing-template-{}-{}",
            std::process::id(),
            NEXT_TEST_PATH.fetch_add(1, Ordering::Relaxed)
        ));
        assert!(matches!(
            CoordinateTemplates::load_templates_from_path(&missing),
            Err(TemplateError::ExternalOpen { .. })
        ));
    }

    #[test]
    fn external_templates_malformed_later_row_does_not_change_existing_registry() {
        let original = CoordinateTemplates::new_default().unwrap();
        let expected_count = original.template_count();
        with_external_rows("C1CC1 |(0,0,;1,0,;0.5,0.866,)|\n[\n", |path| {
            assert!(matches!(
                CoordinateTemplates::load_templates_from_path(path),
                Err(TemplateError::ExternalInvalidSmarts { line: 2, .. })
            ));
        });
        assert_eq!(original.template_count(), expected_count);
    }

    #[test]
    fn registry_templates_set_add_default_preserves_source_prepend_order() {
        let first = "C1CC1 |(0,0,;1,0,;0.5,0.866,)|";
        let second = "N1CC1 |(0,0,;1,0,;0.5,0.866,)|";
        let third = "O1CC1 |(0,0,;1,0,;0.5,0.866,)|";
        let square = "C1CCC1 |(0,0,;1,0,;1,1,;0,1,)|";
        let mut registry = CoordinateTemplates::new_default().unwrap();
        assert_eq!(registry.template_count(), 578);
        with_external_rows(&format!("{first}\n"), |path| {
            registry.set_ring_system_templates(path).unwrap();
        });
        assert_eq!(registry.template_count(), 1);
        assert!(registry.has_template_of_size(3));
        assert!(!registry.has_template_of_size(4));
        with_external_rows(&format!("{second}\n{square}\n{third}\n"), |path| {
            registry.add_ring_system_templates(path).unwrap();
        });
        assert_eq!(registry.template_count(), 4);
        let expected = [second, third, first];
        let actual = registry.templates_of_size(3).unwrap();
        for (row, template) in expected.into_iter().zip(actual) {
            assert_eq!(
                template.query,
                parse_smarts(row, &SmartsParseParams::default()).unwrap()
            );
        }
        assert_eq!(registry.templates_of_size(4).unwrap().len(), 1);
        registry.load_default_templates().unwrap();
        let expected_default = CoordinateTemplates::new_default().unwrap();
        assert_eq!(registry.template_count(), 578);
        assert_eq!(registry.buckets.len(), expected_default.buckets.len());
        for (size, templates) in &registry.buckets {
            let expected = &expected_default.buckets[size];
            assert_eq!(templates.len(), expected.len());
            for (actual, reference) in templates.iter().zip(expected) {
                assert_eq!(actual.query, reference.query);
            }
        }
    }

    #[test]
    fn registry_templates_missing_lookup_and_failed_transactions_follow_source() {
        let first = "C1CC1 |(0,0,;1,0,;0.5,0.866,)|";
        let mut registry = CoordinateTemplates::default();
        assert!(!registry.has_template_of_size(99));
        assert!(registry.matching_templates(99).is_empty());
        assert!(registry.has_template_of_size(99));
        with_external_rows(&format!("{first}\n"), |path| {
            registry.set_ring_system_templates(path).unwrap();
        });
        assert!(!registry.has_template_of_size(99));
        let before = registry.templates_of_size(3).unwrap()[0].query.clone();
        for replacement in [false, true] {
            with_external_rows(&format!("{first}\n[\n"), |path| {
                let result = if replacement {
                    registry.set_ring_system_templates(path)
                } else {
                    registry.add_ring_system_templates(path)
                };
                assert!(matches!(
                    result,
                    Err(TemplateError::ExternalInvalidSmarts { line: 2, .. })
                ));
            });
            assert_eq!(registry.template_count(), 1);
            assert_eq!(registry.templates_of_size(3).unwrap()[0].query, before);
        }
        with_external_rows("", |path| registry.add_ring_system_templates(path).unwrap());
        assert_eq!(registry.template_count(), 1);
        with_external_rows("", |path| registry.set_ring_system_templates(path).unwrap());
        assert_eq!(registry.template_count(), 0);
    }
}
