//! Pinned RDKit ring-system coordinate templates over canonical query graphs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cosmolkit_core::{
    PathError, RingFindingError, RingInfo, connected_components,
    symmetrize_sssr_with_options_from_parts,
};
use cosmolkit_model::{
    QueryAtom, QueryAtomConversionError, QueryGraph, TopologyBlock, TopologyValidationError,
};
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
    NonElementIdentity {
        index: usize,
        source: QueryAtomConversionError,
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
            .map(QueryAtom::try_to_atom)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|source| TemplateError::NonElementIdentity { index, source })?;
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

    #[test]
    fn recovery_dep03_source_asset_rows_order_and_changed_query_predicates() {
        let rows = DEFAULT_TEMPLATE_ROWS.lines().collect::<Vec<_>>();
        assert_eq!(rows.len(), 578);
        assert_eq!(
            rows[0],
            r########"[!#200&D2]1~[!#200&D2]~[!#200]2~[!#200]~[!#200]~[!#200]~1~[!#200]1~[!#200]~[!#200]~[!#200](~[!#200&D2]~[!#200&D2]~1)~[!#200]1~[!#200]~[!#200]~[!#200](~[!#200&D2]~[!#200&D2]~1)~[!#200]1~[!#200&D2]~[!#200&D2]~[!#200](~[!#200]~[!#200]~1)~[!#200]1~[!#200]~[!#200]~[!#200](~[!#200&D2]~[!#200&D2]~1)~[!#200]1~[!#200&D2]~[!#200&D2]~[!#200]~2~[!#200]~[!#200]~1 |(2.33,1.43,;3.56,2.14,;3.56,3.56,;2.33,4.28,;1.1,3.57,;1.1,2.14,;0.12,0.44,;-1.12,-0.27,;-1.12,-1.69,;0.12,-2.4,;1.35,-1.69,;1.35,-0.27,;1.1,-4.1,;1.1,-5.53,;2.33,-6.24,;3.56,-5.52,;3.56,-4.1,;2.33,-3.39,;5.53,-5.52,;5.53,-4.1,;6.76,-3.39,;8,-4.11,;7.99,-5.53,;6.76,-6.24,;8.98,-2.4,;10.21,-1.69,;10.21,-0.27,;8.98,0.44,;7.75,-0.27,;7.75,-1.69,;8,2.14,;6.76,1.43,;5.53,2.14,;5.53,3.56,;6.76,4.27,;7.99,3.56,)|"########
        );
        assert_eq!(
            rows[109],
            r########"[!#200]1~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~1 |(3.14481,-3.43775,;1.52309,-3.87284,;0.132265,-2.96287,;-1.25875,-3.87255,;-2.88037,-3.43713,;-3.59659,-1.97534,;-2.87133,-0.600979,;-3.74061,0.746997,;-3.02553,2.15605,;-1.45178,2.2201,;-0.66539,3.59093,;0.931235,3.59078,;1.71736,2.2198,;3.29112,2.15544,;4.00592,0.746236,;3.13636,-0.601575,;3.86135,-1.97609,)|"########
        );
        assert_eq!(
            rows[110],
            r########"[!#200]1~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~1 |(-1.25916,-3.49323,;0.13675,-2.56295,;1.53542,-3.48915,;3.1713,-3.05689,;3.89741,-1.57913,;3.17124,-0.172643,;4.09547,1.17221,;3.43683,2.64493,;1.83978,2.81498,;0.942355,1.50367,;-0.684476,1.50142,;-1.58566,2.81022,;-3.18206,2.63526,;-3.83568,1.16037,;-2.9066,-0.18134,;-3.62773,-1.59048,;-2.89652,-3.06598,)|"########
        );
        assert_eq!(
            rows[159],
            r########"[!#200&D2]1~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~1 |(3.97333,-1.48791,;3.20544,-2.91896,;1.52915,-3.3031,;0.107207,-2.36163,;-1.31468,-3.30317,;-2.99102,-2.91914,;-3.75905,-1.48815,;-5.41836,-1.51834,;-6.31992,-0.201414,;-5.65569,1.20636,;-4.07651,1.3106,;-3.32641,2.72139,;-1.67355,2.83423,;-0.735358,1.51545,;0.949189,1.5155,;1.88731,2.83432,;3.54018,2.72159,;4.2904,1.31086,;5.86961,1.20677,;6.53401,-0.200923,;5.63262,-1.51794,)|"########
        );
        assert_eq!(
            rows[160],
            r########"[!#200&D2]1~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~1 |(-0.959305,-3.17778,;-1.98477,-4.41288,;-3.61634,-4.05652,;-4.13065,-2.47478,;-3.01012,-1.32982,;-3.20316,0.324015,;-4.55098,1.14682,;-4.48971,2.7944,;-3.06887,3.65712,;-1.6159,2.9395,;-0.199836,3.78552,;1.236,2.9735,;2.67144,3.72556,;4.11246,2.89698,;4.21302,1.25133,;2.88522,0.396593,;2.73171,-1.26138,;3.87925,-2.37927,;3.40283,-3.97282,;1.78024,-4.368,;0.725609,-3.1577,)|"########
        );
        assert_eq!(
            rows[161],
            r########"[!#200&D2]1~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~1 |(0.106988,-2.65765,;1.53681,-3.59266,;3.21588,-3.1907,;3.97423,-1.74386,;5.6427,-1.74849,;6.54193,-0.411921,;5.86808,0.990482,;4.30693,1.0487,;3.47718,2.40474,;1.84356,2.33492,;0.939774,3.66905,;-0.723781,3.66918,;-1.6278,2.33519,;-3.26141,2.40533,;-4.09151,1.04948,;-5.65271,0.991722,;-6.32704,-0.410459,;-5.42828,-1.74728,;-3.75988,-1.74313,;-3.00193,-3.19012,;-1.323,-3.59242,)|"########
        );
        assert_eq!(
            rows[222],
            r########"[!#200]1~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~1 |(6.91103,0.65698,;6.12973,-0.716816,;6.90129,-2.05593,;6.19288,-3.50503,;4.55366,-3.88194,;3.23446,-2.86648,;1.66324,-3.57017,;0.172159,-2.71939,;-1.31998,-3.56816,;-2.88979,-2.86233,;-4.21034,-3.87549,;-5.84797,-3.49498,;-6.55139,-2.04423,;-5.77489,-0.708085,;-6.55121,0.668091,;-5.75296,2.05177,;-4.13239,2.09153,;-3.32037,3.50623,;-1.63942,3.6265,;-0.671284,2.32412,;1.03682,2.32338,;2.00632,3.62493,;3.68711,3.50256,;4.49659,2.08657,;6.11699,2.04302,)|"########
        );
        assert_eq!(
            rows[270],
            r########"[!#200]1~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~1 |(8.01612,10.7082,;6.51203,10.6057,;5.42699,11.7263,;5.5724,13.2123,;4.26,14.0062,;2.9476,13.2123,;3.09301,11.7263,;2.00797,10.6057,;0.503876,10.7082,;-0.31124,9.35638,;0.469584,8.02051,;-0.122302,6.55767,;-1.55508,6.22936,;-2.11219,4.78048,;-1.32989,3.41616,;0.192422,3.17624,;1.17836,4.25753,;2.72861,4.04791,;3.47496,2.71338,;5.04504,2.71338,;5.79139,4.04791,;7.34164,4.25753,;8.32758,3.17624,;9.84989,3.41616,;10.6322,4.78048,;10.0751,6.22936,;8.6423,6.55767,;8.05041,8.02051,;8.83124,9.35638,)|"########
        );
        assert_eq!(
            rows[274],
            r########"[!#200]1~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~1 |(4.6363,-4.7429,;6.32266,-4.38116,;7.14766,-2.94436,;6.50329,-1.52668,;7.4016,-0.182146,;6.65874,1.24663,;7.44376,2.6012,;6.72595,4.08446,;5.07676,4.46147,;3.83492,3.40421,;2.1712,3.88533,;1.02867,2.68067,;-0.716459,2.68084,;-1.85865,3.88576,;-3.52257,3.40525,;-4.7638,4.4632,;-6.41323,4.08745,;-7.13236,2.60497,;-6.34891,1.24972,;-7.09362,-0.177847,;-6.19745,-1.52335,;-6.84394,-2.93978,;-6.02098,-4.37723,;-4.33528,-4.74029,;-2.97763,-3.73934,;-1.3796,-4.47692,;0.150646,-3.63024,;1.68047,-4.47779,;3.2791,-3.74107,)|"########
        );
        assert_eq!(
            rows[277],
            r########"[!#200]1~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~1 |(8.99419,-1.91441,;7.3323,-1.89846,;6.50539,-3.32221,;4.76551,-3.63353,;3.38817,-2.56655,;1.72361,-3.27428,;0.147294,-2.39923,;-1.42884,-3.27417,;-3.09255,-2.56641,;-4.46943,-3.6326,;-6.20747,-3.31951,;-7.0304,-1.89495,;-8.69099,-1.9069,;-9.55729,-0.559464,;-8.83951,0.849329,;-7.23238,0.905406,;-6.34187,2.30078,;-4.61994,2.2491,;-3.64574,3.65196,;-1.84121,3.7468,;-0.76025,2.42733,;1.07498,2.42719,;2.15628,3.74656,;3.96094,3.6513,;4.93456,2.24805,;6.65666,2.29826,;7.54497,0.901664,;9.15167,0.841969,;9.86512,-0.568889,)|"########
        );
        assert_eq!(
            rows[278],
            r########"[!#200]1~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~1 |(1.57759,3.85231,;0.827548,2.55325,;-0.672395,2.55328,;-1.42239,3.85232,;-2.92243,3.8523,;-3.67238,2.55329,;-5.17242,2.55327,;-5.92236,1.25427,;-5.17236,-0.044794,;-5.92241,-1.34383,;-5.17241,-2.64286,;-3.67237,-2.64284,;-2.92237,-3.94188,;-1.42243,-3.9419,;-0.672385,-2.64285,;-1.13589,-1.21626,;0.077649,-0.334596,;1.29116,-1.2163,;0.827559,-2.64287,;1.57756,-3.94191,;3.07756,-3.94191,;3.82754,-2.64288,;5.32757,-2.64286,;6.07757,-1.34386,;5.32757,-0.044794,;6.07757,1.25423,;5.32757,2.55327,;3.82753,2.55325,;3.07753,3.85228,)|"########
        );
        assert_eq!(
            rows[280],
            r########"[!#200]1~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~1 |(8.34997,7.94245,;6.88932,8.18742,;6.11151,9.51452,;6.58804,10.9417,;5.44653,11.9503,;5.55373,13.4257,;4.26,14.1999,;2.96627,13.4257,;3.07347,11.9503,;1.93196,10.9417,;2.40849,9.51452,;1.63068,8.18742,;0.170034,7.94245,;-0.27671,6.46374,;-1.68999,6.09112,;-2.22091,4.66046,;-1.43905,3.32984,;0.062326,3.10605,;1.05954,4.1694,;2.56376,3.85817,;3.49119,5.02174,;5.02881,5.02174,;5.95624,3.85817,;7.46045,4.1694,;8.45767,3.10605,;9.95905,3.32984,;10.7409,4.66046,;10.21,6.09112,;8.79671,6.46374,)|"########
        );
        assert_eq!(
            rows[327],
            r########"[!#200]1~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~1 |(12.2438,1.63048,;11.4172,-0.0981918,;12.2429,-1.82797,;11.3245,-3.54234,;9.32249,-3.78818,;7.96387,-2.41677,;5.81944,-2.35037,;4.39061,-3.72865,;2.12622,-3.87744,;0.220018,-2.75839,;-1.69009,-3.88021,;-3.96413,-3.73377,;-5.40178,-2.35392,;-7.55095,-2.41971,;-8.91338,-3.79012,;-10.9173,-3.54258,;-11.8359,-1.82811,;-11.0102,-0.0982423,;-11.8368,1.63043,;-10.922,3.34313,;-8.95393,3.59341,;-7.63267,2.25073,;-5.57583,2.23199,;-4.25678,3.59555,;-2.16484,3.5761,;-0.846724,2.17802,;1.2536,2.17802,;2.57171,3.5761,;4.66367,3.59555,;5.98274,2.232,;8.03959,2.25074,;9.36086,3.59342,;11.3289,3.34315,)|"########
        );
        assert_eq!(
            rows[336],
            r########"[!#200]1~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~1 |(10.1703,0.416024,;9.42063,-0.97539,;10.1956,-2.32419,;9.45354,-3.77581,;7.77239,-4.10843,;6.50743,-3.00132,;4.75569,-3.4554,;3.57929,-2.18578,;1.7153,-2.01908,;0.20077,-3.04299,;-1.31362,-2.01907,;-3.17713,-2.18591,;-4.3528,-3.45544,;-6.10375,-3.00151,;-7.36833,-4.10808,;-9.04848,-3.77427,;-9.78834,-2.32193,;-9.01115,-0.974207,;-9.75897,0.418179,;-8.92584,1.81415,;-7.26196,1.87576,;-6.45127,3.34371,;-4.70394,3.60287,;-3.53717,2.39435,;-1.71756,2.54341,;-0.715831,3.94878,;1.1287,3.94873,;2.13036,2.54329,;3.95,2.39408,;5.11697,3.60251,;6.86432,3.34311,;7.67462,1.87501,;9.33848,1.8126,)|"########
        );
        assert_eq!(
            rows[577],
            r########"[!#200&D2]1~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200]~[!#200]~[!#200&D2]~[!#200]~[!#200&D2]~[!#200&D2]~[!#200]~[!#200]~[!#200]~[!#200]~[!#200&D2]~1 |(1.22974,-4.97003,;0,-5.67998,;-1.22974,-4.96998,;-2.45952,-5.68002,;-3.68928,-4.97002,;-4.91901,-5.67997,;-6.14876,-4.96997,;-7.37854,-5.68001,;-8.6083,-4.97001,;-9.83802,-5.67996,;-11.0678,-4.96996,;-11.0678,-3.55001,;-12.2976,-2.84001,;-12.2975,-1.41997,;-13.5273,-0.709968,;-13.5273,0.709983,;-12.2975,1.42002,;-11.0678,0.710021,;-9.83806,1.41997,;-8.60831,0.709972,;-7.37853,1.42001,;-7.37856,2.83996,;-6.14878,3.55,;-4.91902,2.84,;-3.68925,3.55004,;-3.68927,4.96999,;-2.4595,5.68003,;-1.22974,4.97003,;0,5.67998,;1.22974,4.96998,;2.45952,5.68002,;2.45949,7.09997,;3.68927,7.81001,;4.91903,7.10001,;4.91901,5.67997,;3.68928,4.97002,;3.68926,3.54998,;4.91901,2.83998,;4.91904,1.42003,;6.1488,0.710028,;6.14877,-0.710011,;4.91905,-1.41996,;4.91902,-2.84,;3.68925,-3.55004,;3.68927,-4.96999,;4.91903,-5.67999,;4.91901,-7.10003,;3.68928,-7.80998,;2.45952,-7.09998,;2.4595,-5.68003,)|"########
        );
        let templates = CoordinateTemplates::new_default().unwrap();
        assert_eq!(templates.template_count(), 578);
    }
    #[test]
    fn recovery_dep03_every_row_has_coordinates_and_preserves_native_query_owner() {
        use cosmolkit_model::{AtomQueryPredicate, BondQueryPredicate, QueryNode};

        // Expectations come from the byte-exact official .6 TemplateSmarts.h
        // rows, not a second invocation of the SMARTS/CX parser under test.
        // The frozen 578-row source uses only [!#200], [!#200&D2/D3/D4/D6],
        // explicit ~ bonds, and ordered (x,y,) coordinates. Reject every
        // unexpected token instead of silently omitting it from expectations.
        let templates = CoordinateTemplates::new_default().expect("all .6 rows compile");
        let rows: Vec<_> = DEFAULT_TEMPLATE_ROWS.lines().collect();
        assert_eq!(rows.len(), 578);
        assert_eq!(templates.template_count(), 578);
        let mut per_size_index = BTreeMap::<usize, usize>::new();
        for (index, row) in rows.into_iter().enumerate() {
            let (graph, coordinates) = row.split_once(" |(").expect("source CX coordinate opener");
            let coordinates = coordinates
                .strip_suffix(")|")
                .expect("source CX coordinate end");
            let mut parts = graph.split('[');
            assert_eq!(parts.next(), Some(""), "source row {index}");
            let mut predicates = Vec::new();
            let mut bond_count = 0;
            for part in parts {
                let (atom, tail) = part.split_once(']').expect("complete source atom bracket");
                assert!(
                    tail.bytes()
                        .all(|byte| byte.is_ascii_digit() || b"~()%".contains(&byte)),
                    "unexpected source graph token in row {index}: {tail}"
                );
                bond_count += tail.bytes().filter(|&byte| byte == b'~').count();
                let not_200 = QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::AtomicNumber(200),
                )));
                let degree = match atom {
                    "!#200" => None,
                    "!#200&D2" => Some(2),
                    "!#200&D3" => Some(3),
                    "!#200&D4" => Some(4),
                    "!#200&D6" => Some(6),
                    other => panic!("unexpected source atom in row {index}: {other}"),
                };
                predicates.push(match degree {
                    None => not_200,
                    Some(value) => QueryNode::And(vec![
                        not_200,
                        QueryNode::Predicate(AtomQueryPredicate::ExplicitDegree(value)),
                    ]),
                });
            }
            assert!(!predicates.is_empty(), "source row {index}");
            let expected_xy: Vec<[f64; 2]> = coordinates
                .split(';')
                .map(|point| {
                    let mut values = point.split(',');
                    let x: f64 = values.next().expect("source x").parse().expect("decimal x");
                    let y: f64 = values.next().expect("source y").parse().expect("decimal y");
                    assert_eq!(values.next(), Some(""), "source empty z in row {index}");
                    assert_eq!(
                        values.next(),
                        None,
                        "extra source coordinate in row {index}"
                    );
                    [x, y]
                })
                .collect();
            assert_eq!(expected_xy.len(), predicates.len(), "source row {index}");

            // The source row's atom sequence is the CX coordinate sequence.
            // Bucket offsets preserve source order within each atom count.
            let slot = per_size_index.entry(predicates.len()).or_default();
            let template = &templates
                .templates_of_size(predicates.len())
                .expect("source bucket")[*slot];
            *slot += 1;
            assert_eq!(template.query.num_atoms(), predicates.len(), "row {index}");
            assert_eq!(template.query.num_bonds(), bond_count, "row {index}");
            for (atom_index, (atom, expected)) in
                template.query.atoms().iter().zip(&predicates).enumerate()
            {
                // Whole-tree equality rejects missing, extra, negated, reordered,
                // OR/XOR or otherwise altered predicate nodes.
                assert_eq!(atom.predicate(), expected, "row {index}, atom {atom_index}");
            }
            for (bond_index, bond) in template.query.bonds().iter().enumerate() {
                assert_eq!(
                    bond.predicate(),
                    &QueryNode::Predicate(BondQueryPredicate::Any),
                    "row {index}, bond {bond_index}"
                );
            }
            // Current query CX lowering stores one source 2D conformer in its
            // 3D coordinate carrier with is_3d=false; inspect that actual owner.
            assert!(template.query.coordinates_2d().is_none(), "row {index}");
            assert_eq!(template.query.conformers_3d().len(), 1, "row {index}");
            let conformer = &template.query.conformers_3d()[0];
            assert!(!conformer.is_3d(), "row {index}");
            assert_eq!(
                conformer.coordinates().len(),
                expected_xy.len(),
                "row {index}"
            );
            for (atom_index, (actual, expected)) in
                conformer.coordinates().iter().zip(&expected_xy).enumerate()
            {
                assert_eq!(
                    actual.map(f64::to_bits),
                    [expected[0], expected[1], 0.0].map(f64::to_bits),
                    "row {index}, coordinate atom {atom_index}"
                );
            }
            assert!(template.rings.is_symm_sssr(), "row {index}");
        }
        assert_eq!(per_size_index.len(), templates.buckets.len());
        for (size, count) in per_size_index {
            assert_eq!(
                templates
                    .templates_of_size(size)
                    .expect("source bucket")
                    .len(),
                count
            );
        }
    }

    #[test]
    fn recovery_dep02_sssr_five_ring_cubane_matches_six_ring_template() {
        let row = "C12C3C4C1C5C2C3C45 |(0,0,;1,0,;1,1,;0,1,;2,2,;3,2,;3,3,;2,3,)|";
        let (template, topology) = CoordinateTemplates::parse_template_row(row, 0).unwrap();
        let expected = if let Some(xy) = template.query.coordinates_2d() {
            xy.to_vec()
        } else {
            template.query.conformers_3d()[0]
                .coordinates()
                .iter()
                .map(|p| [p[0], p[1]])
                .collect()
        };
        let rings = cosmolkit_core::find_sssr(&topology, &Default::default()).unwrap();
        assert_eq!(rings.num_rings(), 5);
        assert_eq!(template.rings.num_rings(), 6);
        let fused = rings
            .atom_rings()
            .iter()
            .map(|r| r.iter().map(|a| a.index()).collect())
            .collect::<Vec<Vec<usize>>>();
        let mut templates = CoordinateTemplates::default();
        templates.buckets.insert(8, vec![template]);
        let frag = crate::embedded_frag::EmbeddedFrag::from_fused_rings(
            &topology,
            &rings,
            &fused,
            true,
            &mut templates,
        )
        .unwrap();
        assert_eq!(frag.atoms.len(), 8);
        assert!(frag.atoms.values().all(|a| a.fixed));
        let mut actual = frag.atoms.values().map(|a| a.loc).collect::<Vec<_>>();
        let mut expected = expected;
        actual.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        expected.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        assert_eq!(actual, expected);
    }
    #[test]
    fn recovery_dep02_same_atom_bond_counts_still_reject_degree_mismatch() {
        // Private registry seam: the path/star are not valid natural ring templates.
        let (_, topology) =
            CoordinateTemplates::parse_template_row("CCCC |(0,0,;1,0,;2,0,;3,0,)|", 0).unwrap();
        let (template, _) =
            CoordinateTemplates::parse_template_row("CC(C)C |(0,0,;1,0,;1,1,;2,0,)|", 1).unwrap();
        let rings = cosmolkit_core::find_sssr(&topology, &Default::default()).unwrap();
        let mut templates = CoordinateTemplates::default();
        templates.buckets.insert(4, vec![template]);
        let mut frag =
            crate::embedded_frag::EmbeddedFrag::from_single(0, &topology, &rings).unwrap();
        let before = format!("{:?}", frag.atoms);
        assert!(
            !frag
                .match_to_template(&[0, 1, 2, 3], &mut templates)
                .unwrap()
        );
        assert_eq!(format!("{:?}", frag.atoms), before);
    }

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
    fn non_element_template_identity_preserves_typed_row_error() {
        match CoordinateTemplates::parse_template_row("[#119]", 7) {
            Err(TemplateError::NonElementIdentity { index, source }) => {
                assert_eq!(index, 7);
                assert_eq!(
                    source,
                    QueryAtomConversionError::NonElementAtomicNumber {
                        atom: cosmolkit_model::AtomId::new(0),
                        atomic_number: 119,
                    }
                );
            }
            Err(error) => panic!("expected typed non-Element identity error, got {error:?}"),
            Ok(_) => panic!("non-Element query carrier unexpectedly became a topology"),
        };

        assert!(
            CoordinateTemplates::parse_template_row("C1CC1 |(0,0,;1,0,;0.5,0.866,)|", 6).is_ok(),
            "ordinary Element-backed template remains representable"
        );
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
