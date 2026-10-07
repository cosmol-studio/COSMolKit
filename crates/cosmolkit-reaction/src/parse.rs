use crate::{Reaction, ReactionParseError, ReactionRole};
use cosmolkit_model::{
    AtomQueryPredicate, BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryNode,
    replace_query_substance_groups,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionParseParams {
    pub use_smiles: bool,
    pub sanitize: bool,
    pub replacements: BTreeMap<String, String>,
    pub allow_cxsmiles: bool,
    /// Declared by RDKit but not consumed by parseReaction. This flag does
    /// not suppress reaction-global CX exceptions.
    pub strict_cxsmiles: bool,
}
impl Default for ReactionParseParams {
    fn default() -> Self {
        // RDKit❗✔️:   bool sanitize = false; /**< sanitize the molecules after building them */
        // RDKit❗✔️:   bool allowCXSMILES = true; /**< recognize and parse CXSMILES*/
        // RDKit❗✔️:   bool strictCXSMILES =
        // RDKit❗✔️:       true; /**< throw an exception if the CXSMILES parsing fails */
        // use_smiles selects the two source overloads; default is SMARTS.
        Self {
            use_smiles: false,
            sanitize: false,
            replacements: BTreeMap::new(),
            allow_cxsmiles: true,
            strict_cxsmiles: true,
        }
    }
}

fn trim_ascii(text: &str) -> &str {
    text.trim_matches(|c: char| c.is_ascii_whitespace())
}
fn slice(text: &str, start: usize, end: usize) -> Result<&str, ReactionParseError> {
    text.get(start..end)
        .ok_or(ReactionParseError::ComponentBounds { start, end })
}
fn separator_positions(text: &str) -> Vec<usize> {
    text.as_bytes()
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| {
            (c == b'>' && (i == 0 || text.as_bytes()[i - 1] != b'-')).then_some(i)
        })
        .collect()
}

fn remove_spaces_around(text: &mut String, pos: usize) {
    // RDKit❗✔️: void removeSpacesAround(std::string &text, size_t pos) {
    // RDKit❗✔️:   auto nextp = pos + 1;
    // RDKit❗✔️:   while (nextp < text.size() && (text[nextp] == ' ' || text[nextp] == '\t')) {
    // RDKit❗✔️:     text.erase(nextp, 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (pos > 0) {
    // RDKit❗✔️:     nextp = pos - 1;
    // RDKit❗✔️:     while (text[nextp] == ' ' || text[nextp] == '\t') {
    // RDKit❗✔️:       text.erase(nextp, 1);
    // RDKit❗✔️:       if (nextp > 0) {
    // RDKit❗✔️:         --nextp;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:
    let next = pos + 1;
    while next < text.len() && matches!(text.as_bytes()[next], b' ' | b'\t') {
        text.remove(next);
    }
    if pos > 0 {
        let mut next = pos - 1;
        while matches!(text.as_bytes()[next], b' ' | b'\t') {
            text.remove(next);
            if next > 0 {
                next -= 1;
            } else {
                break;
            }
        }
    }
}

fn split_components(text: &str) -> Result<Vec<&str>, ReactionParseError> {
    // RDKit❗✔️: std::vector<std::string> splitSmartsIntoComponents(
    // RDKit❗✔️:     const std::string &reactText) {
    // RDKit❗✔️:   std::vector<std::string> res;
    // RDKit❗✔️:   unsigned int pos = 0;
    // RDKit❗✔️:   unsigned int blockStart = 0;
    // RDKit❗✔️:   unsigned int level = 0;
    // RDKit❗✔️:   unsigned int inBlock = 0;
    // RDKit❗✔️:   while (pos < reactText.size()) {
    // RDKit❗✔️:     if (reactText[pos] == '(') {
    // RDKit❗✔️:       if (pos == blockStart) {
    // RDKit❗✔️:         inBlock = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       ++level;
    // RDKit❗✔️:     } else if (reactText[pos] == ')') {
    // RDKit❗✔️:       if (level == 1 && inBlock) {
    // RDKit❗✔️:         // this closes a block
    // RDKit❗✔️:         inBlock = 2;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       --level;
    // RDKit❗✔️:     } else if (level == 0 && reactText[pos] == '.') {
    // RDKit❗✔️:       if (inBlock == 2) {
    // RDKit❗✔️:         std::string element =
    // RDKit❗✔️:             reactText.substr(blockStart + 1, pos - blockStart - 2);
    // RDKit❗✔️:         res.push_back(element);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         std::string element = reactText.substr(blockStart, pos - blockStart);
    // RDKit❗✔️:         res.push_back(element);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       blockStart = pos + 1;
    // RDKit❗✔️:       inBlock = 0;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++pos;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (blockStart < pos) {
    // RDKit❗✔️:     if (inBlock == 2) {
    // RDKit❗✔️:       std::string element =
    // RDKit❗✔️:           reactText.substr(blockStart + 1, pos - blockStart - 2);
    // RDKit❗✔️:       res.push_back(element);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       std::string element = reactText.substr(blockStart, pos - blockStart);
    // RDKit❗✔️:       res.push_back(element);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    if text.len() > u32::MAX as usize {
        return Err(ReactionParseError::OffsetOverflow);
    }
    let mut result = Vec::new();
    let mut start = 0;
    let mut level = 0u32;
    let mut block = 0;
    for (pos, &byte) in text.as_bytes().iter().enumerate() {
        match byte {
            b'(' => {
                if pos == start {
                    block = 1;
                }
                level = level.wrapping_add(1);
            }
            b')' => {
                if level == 1 && block != 0 {
                    block = 2;
                }
                level = level.wrapping_sub(1);
            }
            b'.' if level == 0 => {
                result.push(if block == 2 {
                    slice(text, start + 1, pos.saturating_sub(1))?
                } else {
                    slice(text, start, pos)?
                });
                start = pos + 1;
                block = 0;
            }
            _ => {}
        }
    }
    if start < text.len() {
        result.push(if block == 2 {
            slice(text, start + 1, text.len().saturating_sub(1))?
        } else {
            slice(text, start, text.len())?
        });
    }
    Ok(result)
}

fn construct_component(
    text: &str,
    params: &ReactionParseParams,
    role: ReactionRole,
    template: usize,
) -> Result<QueryGraph, ReactionParseError> {
    // RDKit❗❌: std::unique_ptr<RWMol> constructMolFromString(
    // RDKit❗❌:     const std::string &txt, const ReactionSmartsParserParams &params,
    // RDKit❗❌:     bool useSmiles) {
    // RDKit❗❌:   if (!useSmiles) {
    // RDKit❗❌:     SmilesParse::SmartsParserParams ps;
    // RDKit❗❌:     ps.replacements = params.replacements;
    // RDKit❗❌:     ps.allowCXSMILES = false;
    // RDKit❗❌:     ps.parseName = false;
    // RDKit❗❌:     ps.mergeHs = false;
    // RDKit❗❌:     ps.skipCleanup = true;
    // RDKit❗❌:     return SmilesParse::MolFromSmarts(txt, ps);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     SmilesParse::SmilesParserParams ps;
    // RDKit❗❌:     ps.replacements = params.replacements;
    // RDKit❗❌:     ps.allowCXSMILES = false;
    // RDKit❗❌:     ps.parseName = false;
    // RDKit❗❌:     ps.sanitize = params.sanitize;
    // RDKit❗❌:     ps.removeHs = false;
    // RDKit❗❌:     ps.skipCleanup = true;
    // RDKit❗❌:     return SmilesParse::MolFromSmiles(txt, ps);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // RDKit❗❌:
    if !params.use_smiles {
        return cosmolkit_search::parse_smarts(
            text,
            &cosmolkit_search::SmartsParseParams {
                replacements: params.replacements.clone(),
                allow_cxsmiles: false,
                parse_name: false,
                merge_hs: false,
                skip_cleanup: true,
                ..Default::default()
            },
        )
        .map_err(|source| ReactionParseError::Smarts {
            role,
            template,
            text: text.into(),
            source,
        });
    }
    let record = cosmolkit_smiles::parse_smiles(
        text,
        &cosmolkit_smiles::SmilesParseParams {
            replacements: params.replacements.clone(),
            allow_cxsmiles: false,
            parse_name: false,
            sanitize: params.sanitize,
            remove_hydrogens: false,
            skip_cleanup: true,
            ..Default::default()
        },
    )
    .map_err(|source| ReactionParseError::Smiles {
        role,
        template,
        text: text.into(),
        source,
    })?;
    // Canonical carrier projection preserves hasQuery=false and all detached
    // properties. No live Molecule or alternative query AST is constructed.
    let atoms = record
        .topology
        .atoms
        .into_iter()
        .map(|atom| {
            let predicate = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(
                atom.element().atomic_number(),
            ));
            QueryAtom::from_carrier_parts(atom, predicate)
        })
        .collect();
    let bonds = record
        .topology
        .bonds
        .into_iter()
        .map(|bond| {
            let predicate = QueryNode::predicate(BondQueryPredicate::Order(bond.order()));
            QueryBond::from_carrier_parts(bond, predicate)
        })
        .collect();
    let mut graph = QueryGraph::from_parts(
        atoms,
        bonds,
        record.properties.props().clone(),
        record.coordinates.conformers_2d,
        record.coordinates.conformers_3d,
        record.topology.stereo_groups,
    )
    .map_err(|source| ReactionParseError::Model {
        role,
        template,
        source,
    })?;
    replace_query_substance_groups(&mut graph, record.topology.substance_groups).map_err(
        |source| ReactionParseError::Model {
            role,
            template,
            source,
        },
    )?;
    Ok(graph)
}

pub fn parse_smirks_with_params(
    orig_text: &str,
    params: &ReactionParseParams,
) -> Result<Reaction, ReactionParseError> {
    // RDKit❗❌: std::unique_ptr<ChemicalReaction> parseReaction(
    // RDKit❗❌:     const std::string &origText, const ReactionSmartsParserParams &params,
    // RDKit❗❌:     bool useSmiles) {
    // RDKit❗❌:   std::string text = origText;
    // RDKit❗❌:   std::string cxPart;
    // RDKit❗❌:   if (params.allowCXSMILES) {
    // RDKit❗❌:     auto sidx = origText.find_first_of("|");
    // RDKit❗❌:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit❗❌:       text = origText.substr(0, sidx);
    // RDKit❗❌:       cxPart = boost::trim_copy(origText.substr(sidx, origText.size() - sidx));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // remove any spaces at the beginning, end, or before the '>'s
    // RDKit❗❌:   boost::trim(text);
    // RDKit❗❌:   std::vector<std::size_t> pos;
    // RDKit❗❌:   for (std::size_t i = 0; i < text.length(); ++i) {
    // RDKit❗❌:     if (text[i] == '>' && (i == 0 || text[i - 1] != '-')) {
    // RDKit❗❌:       pos.push_back(i);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (pos.size() < 2) {
    // RDKit❗❌:     throw ChemicalReactionParserException(
    // RDKit❗❌:         "a reaction requires at least two > characters");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // remove spaces around ">" symbols
    // RDKit❗❌:   for (auto p : boost::make_iterator_range(pos.rbegin(), pos.rend())) {
    // RDKit❗❌:     removeSpacesAround(text, p);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // remove spaces around "." symbols
    // RDKit❗❌:   pos.clear();
    // RDKit❗❌:   for (std::size_t i = 0; i < text.length(); ++i) {
    // RDKit❗❌:     if (text[i] == '.') {
    // RDKit❗❌:       pos.push_back(i);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto p : boost::make_iterator_range(pos.rbegin(), pos.rend())) {
    // RDKit❗❌:     removeSpacesAround(text, p);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // we shouldn't have whitespace left in the reaction string, so go ahead and
    // RDKit❗❌:   // split and strip:
    // RDKit❗❌:   auto sidx = text.find_first_of(" \t");
    // RDKit❗❌:   if (sidx != std::string::npos && sidx != 0) {
    // RDKit❗❌:     text = text.substr(0, sidx);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // re-find the '>' characters so that we can split on them
    // RDKit❗❌:   pos.clear();
    // RDKit❗❌:   for (std::size_t i = 0; i < text.length(); ++i) {
    // RDKit❗❌:     if (text[i] == '>' && (i == 0 || text[i - 1] != '-')) {
    // RDKit❗❌:       pos.push_back(i);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // there's always the chance that one or more of the ">" was in the name
    // RDKit❗❌:   // part, so verify that we have exactly two:
    // RDKit❗❌:   if (pos.size() < 2) {
    // RDKit❗❌:     throw ChemicalReactionParserException(
    // RDKit❗❌:         "a reaction requires at least two > characters");
    // RDKit❗❌:   }
    // RDKit❗❌:   if (pos.size() > 2) {
    // RDKit❗❌:     throw ChemicalReactionParserException("multi-step reactions not supported");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto pos1 = pos[0];
    // RDKit❗❌:   auto pos2 = pos[1];
    // RDKit❗❌:
    // RDKit❗❌:   auto reactText = text.substr(0, pos1);
    // RDKit❗❌:   std::string agentText;
    // RDKit❗❌:   if (pos2 != pos1 + 1) {
    // RDKit❗❌:     agentText = text.substr(pos1 + 1, (pos2 - pos1) - 1);
    // RDKit❗❌:   }
    // RDKit❗❌:   auto productText = text.substr(pos2 + 1);
    // RDKit❗❌:
    // RDKit❗❌:   // recognize changes within the same molecules, e.g., intra molecular bond
    // RDKit❗❌:   // formation therefore we need to correctly interpret parenthesis and dots
    // RDKit❗❌:   // in the reaction smarts
    // RDKit❗❌:   auto reactSmarts = DaylightParserUtils::splitSmartsIntoComponents(reactText);
    // RDKit❗❌:   auto productSmarts =
    // RDKit❗❌:       DaylightParserUtils::splitSmartsIntoComponents(productText);
    // RDKit❗❌:
    // RDKit❗❌:   auto rxn = std::make_unique<ChemicalReaction>();
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &txt : reactSmarts) {
    // RDKit❗❌:     auto mol =
    // RDKit❗❌:         DaylightParserUtils::constructMolFromString(txt, params, useSmiles);
    // RDKit❗❌:     if (!mol) {
    // RDKit❗❌:       std::string errMsg = "Problems constructing reactant from SMARTS: ";
    // RDKit❗❌:       errMsg += txt;
    // RDKit❗❌:       throw ChemicalReactionParserException(errMsg);
    // RDKit❗❌:     }
    // RDKit❗❌:     rxn->addReactantTemplate(ROMOL_SPTR(mol.release()));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &txt : productSmarts) {
    // RDKit❗❌:     auto mol =
    // RDKit❗❌:         DaylightParserUtils::constructMolFromString(txt, params, useSmiles);
    // RDKit❗❌:     if (!mol) {
    // RDKit❗❌:       std::string errMsg = "Problems constructing product from SMARTS: ";
    // RDKit❗❌:       errMsg += txt;
    // RDKit❗❌:       throw ChemicalReactionParserException(errMsg);
    // RDKit❗❌:     }
    // RDKit❗❌:     rxn->addProductTemplate(ROMOL_SPTR(mol.release()));
    // RDKit❗❌:   }
    // RDKit❗❌:   updateProductsStereochem(rxn.get());
    // RDKit❗❌:
    // RDKit❗❌:   // allow a reaction template to have no agent specified
    // RDKit❗❌:   if (agentText.size() != 0) {
    // RDKit❗❌:     auto agentMol = DaylightParserUtils::constructMolFromString(
    // RDKit❗❌:         agentText, params, useSmiles);
    // RDKit❗❌:     if (!agentMol) {
    // RDKit❗❌:       std::string errMsg = "Problems constructing agent from SMARTS: ";
    // RDKit❗❌:       errMsg += agentText;
    // RDKit❗❌:       throw ChemicalReactionParserException(errMsg);
    // RDKit❗❌:     }
    // RDKit❗❌:     std::vector<ROMOL_SPTR> agents = MolOps::getMolFrags(*agentMol, false);
    // RDKit❗❌:     for (auto &agent : agents) {
    // RDKit❗❌:       rxn->addAgentTemplate(agent);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (params.allowCXSMILES && !cxPart.empty()) {
    // RDKit❗❌:     unsigned int startAtomIdx = 0;
    // RDKit❗❌:     unsigned int startBondIdx = 0;
    // RDKit❗❌:     for (auto &mol : boost::make_iterator_range(rxn->beginReactantTemplates(),
    // RDKit❗❌:                                                 rxn->endReactantTemplates())) {
    // RDKit❗❌:       SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
    // RDKit❗❌:                                         cxPart, startAtomIdx, startBondIdx);
    // RDKit❗❌:       startAtomIdx += mol->getNumAtoms();
    // RDKit❗❌:       startBondIdx += mol->getNumBonds();
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto &mol : boost::make_iterator_range(rxn->beginAgentTemplates(),
    // RDKit❗❌:                                                 rxn->endAgentTemplates())) {
    // RDKit❗❌:       SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
    // RDKit❗❌:                                         cxPart, startAtomIdx, startBondIdx);
    // RDKit❗❌:       startAtomIdx += mol->getNumAtoms();
    // RDKit❗❌:       startBondIdx += mol->getNumBonds();
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto &mol : boost::make_iterator_range(rxn->beginProductTemplates(),
    // RDKit❗❌:                                                 rxn->endProductTemplates())) {
    // RDKit❗❌:       SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
    // RDKit❗❌:                                         cxPart, startAtomIdx, startBondIdx);
    // RDKit❗❌:       startAtomIdx += mol->getNumAtoms();
    // RDKit❗❌:       startBondIdx += mol->getNumBonds();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // final cleanups:
    // RDKit❗❌:   for (auto &mol : boost::make_iterator_range(rxn->beginReactantTemplates(),
    // RDKit❗❌:                                               rxn->endReactantTemplates())) {
    // RDKit❗❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto &mol : boost::make_iterator_range(rxn->beginAgentTemplates(),
    // RDKit❗❌:                                               rxn->endAgentTemplates())) {
    // RDKit❗❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto &mol : boost::make_iterator_range(rxn->beginProductTemplates(),
    // RDKit❗❌:                                               rxn->endProductTemplates())) {
    // RDKit❗❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // "SMARTS"-based reactions have implicit properties
    // RDKit❗❌:   rxn->setImplicitPropertiesFlag(true);
    // RDKit❗❌:
    // RDKit❗❌:   return rxn;
    // RDKit❗❌: }
    let (text, cx) = if params.allow_cxsmiles {
        if let Some(pos) = orig_text.find('|').filter(|&pos| pos != 0) {
            (&orig_text[..pos], trim_ascii(&orig_text[pos..]))
        } else {
            (orig_text, "")
        }
    } else {
        (orig_text, "")
    };
    let mut text = trim_ascii(text).to_owned();
    let positions = separator_positions(&text);
    if positions.len() < 2 {
        return Err(ReactionParseError::Separators {
            count: positions.len(),
        });
    }
    for pos in positions.into_iter().rev() {
        remove_spaces_around(&mut text, pos);
    }
    let dots: Vec<_> = text
        .as_bytes()
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| (c == b'.').then_some(i))
        .collect();
    for pos in dots.into_iter().rev() {
        remove_spaces_around(&mut text, pos);
    }
    if let Some(pos) = text.find([' ', '\t']).filter(|&pos| pos != 0) {
        text.truncate(pos);
    }
    let positions = separator_positions(&text);
    if positions.len() < 2 {
        return Err(ReactionParseError::Separators {
            count: positions.len(),
        });
    }
    if positions.len() > 2 {
        return Err(ReactionParseError::MultiStep {
            count: positions.len(),
        });
    }
    let (first, second) = (positions[0], positions[1]);
    let mut reaction = Reaction::new();
    for (template, text) in split_components(&text[..first])?.into_iter().enumerate() {
        reaction.reactants.push(construct_component(
            text,
            params,
            ReactionRole::Reactant,
            template,
        )?);
    }
    for (template, text) in split_components(&text[second + 1..])?
        .into_iter()
        .enumerate()
    {
        reaction.products.push(construct_component(
            text,
            params,
            ReactionRole::Product,
            template,
        )?);
    }
    crate::template_stereo::update_products_stereochem(&mut reaction)?;
    let agent_text = &text[first + 1..second];
    if !agent_text.is_empty() {
        let graph = construct_component(agent_text, params, ReactionRole::Agent, 0)?;
        reaction.agents = cosmolkit_search::query_graph_fragments(&graph)
            .map_err(|source| ReactionParseError::AgentFragments { source })?;
    }
    let mut start_atom = 0usize;
    let mut start_bond = 0usize;
    if !cx.is_empty() {
        for (role, templates) in [
            (ReactionRole::Reactant, &mut reaction.reactants),
            (ReactionRole::Agent, &mut reaction.agents),
            (ReactionRole::Product, &mut reaction.products),
        ] {
            for (template, graph) in templates.iter_mut().enumerate() {
                let parsed = cosmolkit_cx::parse_cx_extensions_with_atom_window(
                    cx,
                    start_atom,
                    graph.num_atoms(),
                )
                .map_err(|source| ReactionParseError::CxParse {
                    role,
                    template,
                    start_atom,
                    start_bond,
                    source,
                })?;
                cosmolkit_search::apply_cx_to_query_graph_with_offsets(
                    graph, &parsed, start_atom, start_bond,
                )
                .map_err(|source| ReactionParseError::CxLowering {
                    role,
                    template,
                    start_atom,
                    start_bond,
                    source,
                })?;
                start_atom = start_atom
                    .checked_add(graph.num_atoms())
                    .filter(|&n| n <= u32::MAX as usize)
                    .ok_or(ReactionParseError::OffsetOverflow)?;
                start_bond = start_bond
                    .checked_add(graph.num_bonds())
                    .filter(|&n| n <= u32::MAX as usize)
                    .ok_or(ReactionParseError::OffsetOverflow)?;
            }
        }
    }
    // RDKit✔️❌:   // final cleanups:
    // RDKit✔️❌:   for (auto &mol : boost::make_iterator_range(rxn->beginReactantTemplates(),
    // RDKit✔️❌:                                               rxn->endReactantTemplates())) {
    // RDKit✔️❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto &mol : boost::make_iterator_range(rxn->beginAgentTemplates(),
    // RDKit✔️❌:                                               rxn->endAgentTemplates())) {
    // RDKit✔️❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto &mol : boost::make_iterator_range(rxn->beginProductTemplates(),
    // RDKit✔️❌:                                               rxn->endProductTemplates())) {
    // RDKit✔️❌:     SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
    // RDKit✔️❌:   }
    // Keep native reactant, agent, product cleanup order and stop at the
    // first reached structural error, retaining its typed source and location.
    for (role, graphs) in [
        (ReactionRole::Reactant, &mut reaction.reactants),
        (ReactionRole::Agent, &mut reaction.agents),
        (ReactionRole::Product, &mut reaction.products),
    ] {
        for (template, graph) in graphs.iter_mut().enumerate() {
            cosmolkit_search::cleanup_query_graph_parser_state(graph).map_err(|source| {
                ReactionParseError::ParserCleanup {
                    role,
                    template,
                    source,
                }
            })?;
        }
    }
    reaction.implicit_properties = true;
    Ok(reaction)
}

pub fn parse_smirks(text: &str) -> Result<Reaction, ReactionParseError> {
    // RDKit❗✔️: std::unique_ptr<ChemicalReaction> ReactionFromSmarts(
    // RDKit❗✔️:     const std::string &origText, const ReactionSmartsParserParams &options) {
    // RDKit❗✔️:   return parseReaction(origText, options, false);
    // RDKit❗✔️: }
    parse_smirks_with_params(text, &ReactionParseParams::default())
}
