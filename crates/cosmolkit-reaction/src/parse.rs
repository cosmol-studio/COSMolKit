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

fn trim_source_text(text: &str) -> &str {
    // Reuse the sole counted-byte Boost/C-locale trim owner, including VT.
    // Only ASCII edge bytes are removed, so these boundaries stay valid UTF8.
    let trimmed = cosmolkit_search::trim_source_whitespace(text.as_bytes());
    let start = trimmed.as_ptr() as usize - text.as_ptr() as usize;
    &text[start..start + trimmed.len()]
}
#[cfg(test)]
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

fn remove_spaces_around(text: &mut String, pos: usize) -> Result<(), ReactionParseError> {
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
    let next = pos.wrapping_add(1);
    while next < text.len() && matches!(text.as_bytes()[next], b' ' | b'\t') {
        text.remove(next);
    }
    if pos > 0 {
        let mut next = pos - 1;
        while matches!(source_string_byte(text, next)?, b' ' | b'\t') {
            text.remove(next);
            if next > 0 {
                next -= 1;
            } else {
                break;
            }
        }
    }
    Ok(())
}

fn source_string_byte(text: &str, pos: usize) -> Result<u8, ReactionParseError> {
    // libstdc++❗✔️:       _GLIBCXX_NODISCARD _GLIBCXX20_CONSTEXPR
    // libstdc++❗✔️:       reference
    // libstdc++❗✔️:       operator[](size_type __pos)
    // libstdc++❗✔️:       {
    // libstdc++❗✔️:         // Allow pos == size() both in C++98 mode, as v3 extension,
    // libstdc++❗✔️: 	// and in C++11 mode.
    // libstdc++❗✔️: 	__glibcxx_assert(__pos <= size());
    // libstdc++❗✔️:         // In pedantic mode be strict in C++98 mode.
    // libstdc++❗✔️: 	_GLIBCXX_DEBUG_PEDASSERT(__cplusplus >= 201103L || __pos < size());
    // libstdc++❗✔️: 	return _M_data()[__pos];
    // libstdc++❗✔️:       }
    // libstdc++❗✔️:
    // libstdc++❗✔️:       /**
    // basic_string stores a NUL at size(); that defined endpoint is not an
    // arbitrary out-of-range fallback. Larger indexes violate its assertion.
    if pos == text.len() {
        Ok(0)
    } else {
        text.as_bytes()
            .get(pos)
            .copied()
            .ok_or(ReactionParseError::ComponentBounds {
                start: pos,
                end: pos,
            })
    }
}

fn component_substr_source(text: &[u8], pos: u32, count: u32) -> Result<&[u8], ReactionParseError> {
    // libstdc++❗🔝:       _GLIBCXX20_CONSTEXPR
    // libstdc++❗🔝:       size_type
    // libstdc++❗🔝:       _M_check(size_type __pos, const char* __s) const
    // libstdc++❗🔝:       {
    // libstdc++❗🔝: 	if (__pos > this->size())
    // libstdc++❗🔝: 	  __throw_out_of_range_fmt(__N("%s: __pos (which is %zu) > "
    // libstdc++❗🔝: 				       "this->size() (which is %zu)"),
    // libstdc++❗🔝: 				   __s, __pos, this->size());
    // libstdc++❗🔝: 	return __pos;
    // libstdc++❗🔝:       }
    // libstdc++❗🔝:       _GLIBCXX20_CONSTEXPR
    // libstdc++❗🔝:       size_type
    // libstdc++❗🔝:       _M_limit(size_type __pos, size_type __off) const _GLIBCXX_NOEXCEPT
    // libstdc++❗🔝:       {
    // libstdc++❗🔝: 	const bool __testoff =  __off < this->size() - __pos;
    // libstdc++❗🔝: 	return __testoff ? __off : this->size() - __pos;
    // libstdc++❗🔝:       }
    // libstdc++❗🔝:       _GLIBCXX20_CONSTEXPR
    // libstdc++❗🔝:       basic_string(const basic_string& __str, size_type __pos,
    // libstdc++❗🔝: 		   size_type __n)
    // libstdc++❗🔝:       : _M_dataplus(_M_local_data())
    // libstdc++❗🔝:       {
    // libstdc++❗🔝: 	const _CharT* __start = __str._M_data()
    // libstdc++❗🔝: 	  + __str._M_check(__pos, "basic_string::basic_string");
    // libstdc++❗🔝: 	_M_construct(__start, __start + __str._M_limit(__pos, __n),
    // libstdc++❗🔝: 		     std::forward_iterator_tag());
    // libstdc++❗🔝:       }
    // libstdc++❗🔝:       _GLIBCXX_NODISCARD _GLIBCXX20_CONSTEXPR
    // libstdc++❗🔝:       basic_string
    // libstdc++❗🔝:       substr(size_type __pos = 0, size_type __n = npos) const
    // libstdc++❗🔝:       { return basic_string(*this,
    // libstdc++❗🔝: 			    _M_check(__pos, "basic_string::substr"), __n); }
    // Native substr checks only the start, then limits count to the remaining
    // counted bytes. A borrowed byte view preserves those bytes, including NUL
    // and non-UTF8, without the native substring allocation/copy: O(1) vs O(L).
    let start = pos as usize;
    if start > text.len() {
        return Err(ReactionParseError::ComponentBounds { start, end: start });
    }
    let remaining = text.len() - start;
    let count = count as usize;
    let length = if count < remaining { count } else { remaining };
    Ok(&text[start..start + length])
}

fn split_components_source(text: &[u8]) -> Result<Vec<&[u8]>, ReactionParseError> {
    // RDKit❗🔝: std::vector<std::string> splitSmartsIntoComponents(
    // RDKit❗🔝:     const std::string &reactText) {
    // RDKit❗🔝:   std::vector<std::string> res;
    // RDKit❗🔝:   unsigned int pos = 0;
    // RDKit❗🔝:   unsigned int blockStart = 0;
    // RDKit❗🔝:   unsigned int level = 0;
    // RDKit❗🔝:   unsigned int inBlock = 0;
    // RDKit❗🔝:   while (pos < reactText.size()) {
    // RDKit❗🔝:     if (reactText[pos] == '(') {
    // RDKit❗🔝:       if (pos == blockStart) {
    // RDKit❗🔝:         inBlock = 1;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       ++level;
    // RDKit❗🔝:     } else if (reactText[pos] == ')') {
    // RDKit❗🔝:       if (level == 1 && inBlock) {
    // RDKit❗🔝:         // this closes a block
    // RDKit❗🔝:         inBlock = 2;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       --level;
    // RDKit❗🔝:     } else if (level == 0 && reactText[pos] == '.') {
    // RDKit❗🔝:       if (inBlock == 2) {
    // RDKit❗🔝:         std::string element =
    // RDKit❗🔝:             reactText.substr(blockStart + 1, pos - blockStart - 2);
    // RDKit❗🔝:         res.push_back(element);
    // RDKit❗🔝:       } else {
    // RDKit❗🔝:         std::string element = reactText.substr(blockStart, pos - blockStart);
    // RDKit❗🔝:         res.push_back(element);
    // RDKit❗🔝:       }
    // RDKit❗🔝:       blockStart = pos + 1;
    // RDKit❗🔝:       inBlock = 0;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     ++pos;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (blockStart < pos) {
    // RDKit❗🔝:     if (inBlock == 2) {
    // RDKit❗🔝:       std::string element =
    // RDKit❗🔝:           reactText.substr(blockStart + 1, pos - blockStart - 2);
    // RDKit❗🔝:       res.push_back(element);
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       std::string element = reactText.substr(blockStart, pos - blockStart);
    // RDKit❗🔝:       res.push_back(element);
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return res;
    // RDKit❗🔝: }
    // Use counted bytes and exact uint32 state/arithmetic. No parenthesis
    // validation, character-boundary repair, trimming or empty-token filtering.
    // One scan and O(K) borrowed result headers replace native O(K+copied bytes)
    // owned substring storage; input lifetime remains borrowed in this private
    // read-only projection. Wider-than-u32 input retains its structural range
    // error because native's unsigned cursor would wrap instead of terminate.
    if text.len() > u32::MAX as usize {
        return Err(ReactionParseError::OffsetOverflow);
    }
    let mut result = Vec::new();
    let mut pos = 0u32;
    let mut block_start = 0u32;
    let mut level = 0u32;
    let mut in_block = 0u32;
    while (pos as usize) < text.len() {
        if text[pos as usize] == b'(' {
            if pos == block_start {
                in_block = 1;
            }
            level = level.wrapping_add(1);
        } else if text[pos as usize] == b')' {
            if level == 1 && in_block != 0 {
                in_block = 2;
            }
            level = level.wrapping_sub(1);
        } else if level == 0 && text[pos as usize] == b'.' {
            result.push(if in_block == 2 {
                component_substr_source(
                    text,
                    block_start.wrapping_add(1),
                    pos.wrapping_sub(block_start).wrapping_sub(2),
                )?
            } else {
                component_substr_source(text, block_start, pos.wrapping_sub(block_start))?
            });
            block_start = pos.wrapping_add(1);
            in_block = 0;
        }
        pos = pos.wrapping_add(1);
    }
    if block_start < pos {
        result.push(if in_block == 2 {
            component_substr_source(
                text,
                block_start.wrapping_add(1),
                pos.wrapping_sub(block_start).wrapping_sub(2),
            )?
        } else {
            component_substr_source(text, block_start, pos.wrapping_sub(block_start))?
        });
    }
    Ok(result)
}

#[cfg(test)]
fn split_components(text: &str) -> Result<Vec<&str>, ReactionParseError> {
    // Existing UTF8 caller projection delegates the sole byte algorithm, then
    // checks exactly the same source byte slices. Never replace malformed
    // output bytes, guess balanced groups or silently discard a component.
    split_components_source(text.as_bytes())?
        .into_iter()
        .map(|component| {
            // Every result borrows this same input allocation, including its end.
            let start = component.as_ptr() as usize - text.as_ptr() as usize;
            slice(text, start, start + component.len())
        })
        .collect()
}

fn construct_component(
    text: impl AsRef<[u8]>,
    params: &ReactionParseParams,
    use_smiles: bool,
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
    // Reuse the sole complete SMARTS/SMILES parser owners with exactly the
    // source overrides. No parser retry, query-to-molecule heuristic, name/CX
    // extraction, H merge/removal or parser cleanup is added here.
    // The SMILES result needs an explicit detached carrier projection in Rust;
    // native returns the same RWMol. Extra graph validation/wrapping/property
    // copies remain a marked cost gap rather than native cost equivalence.
    let text = text.as_ref();
    if !use_smiles {
        return cosmolkit_search::parse_smarts(
            text,
            &cosmolkit_search::SmartsParseParams {
                replacements: params
                    .replacements
                    .iter()
                    .map(|(key, value)| (key.into(), value.into()))
                    .collect(),
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
    // SMARTS consumes source bytes directly. The existing SMILES owner takes
    // str; preserve a non-UTF8 source component and its explicit boundary error
    // rather than replacing/dropping bytes or pretending a native lexer error.
    let smiles_text =
        std::str::from_utf8(text).map_err(|source| ReactionParseError::ComponentEncoding {
            role,
            template,
            text: text.into(),
            source,
        })?;
    let record = cosmolkit_smiles::parse_smiles_complete_source(
        smiles_text,
        &cosmolkit_smiles::SmilesParseParams {
            replacements: params
                .replacements
                .iter()
                .map(|(key, value)| (key.into(), value.into()))
                .collect(),
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
    project_smiles_component_record(record, role, template)
}

fn project_smiles_component_record(
    record: cosmolkit_smiles::SmilesRecord,
    role: ReactionRole,
    template: usize,
) -> Result<QueryGraph, ReactionParseError> {
    // This is the existing detached value projection of the returned RWMol,
    // not another parser/chemical algorithm. Move carrier/coordinate/group
    // values, preserving source hasQuery=false, actual property record order
    // and counted-byte tagged values, and the actual mixed conformer order.
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
        record
            .properties
            .ordered_props()
            .map(|(key, value)| (key.clone(), value.clone())),
        record.coordinates.conformers_2d,
        record.coordinates.conformers_3d,
        record.topology.stereo_groups,
    )
    .map_err(|source| ReactionParseError::Model {
        role,
        template,
        source,
    })?;
    graph
        .set_source_conformer_order(record.coordinates.source_conformer_order)
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
    if params.use_smiles {
        reaction_from_smiles_source(orig_text, params)
    } else {
        reaction_from_smarts_source(orig_text, params)
    }
}

fn parse_reaction_source(
    orig_text: &str,
    params: &ReactionParseParams,
    use_smiles: bool,
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
            (&orig_text[..pos], trim_source_text(&orig_text[pos..]))
        } else {
            (orig_text, "")
        }
    } else {
        (orig_text, "")
    };
    let mut text = trim_source_text(text).to_owned();
    let positions = separator_positions(&text);
    if positions.len() < 2 {
        return Err(ReactionParseError::Separators {
            count: positions.len(),
        });
    }
    for pos in positions.into_iter().rev() {
        remove_spaces_around(&mut text, pos)?;
    }
    let dots: Vec<_> = text
        .as_bytes()
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| (c == b'.').then_some(i))
        .collect();
    for pos in dots.into_iter().rev() {
        remove_spaces_around(&mut text, pos)?;
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
    // Both source component splits complete before any template construction.
    let reactants = split_components_source(text[..first].as_bytes())?;
    let products = split_components_source(text[second + 1..].as_bytes())?;
    let mut reaction = Reaction::new();
    for (template, text) in reactants.into_iter().enumerate() {
        let graph =
            construct_component(text, params, use_smiles, ReactionRole::Reactant, template)?;
        reaction.add_reactant_template_source(graph);
    }
    for (template, text) in products.into_iter().enumerate() {
        let graph = construct_component(text, params, use_smiles, ReactionRole::Product, template)?;
        reaction.add_product_template_source(graph);
    }
    crate::template_stereo::update_products_stereochem(&mut reaction)?;
    let agent_text = &text[first + 1..second];
    if !agent_text.is_empty() {
        let graph = construct_component(agent_text, params, use_smiles, ReactionRole::Agent, 0)?;
        let agents = cosmolkit_search::query_graph_fragments(&graph)
            .map_err(|source| ReactionParseError::AgentFragments { source })?;
        for graph in agents {
            reaction.add_agent_template_source(graph);
        }
    }
    let mut start_atom = 0u32;
    let mut start_bond = 0u32;
    if params.allow_cxsmiles && !cx.is_empty() {
        for (role, templates) in [
            (ReactionRole::Reactant, &mut reaction.reactants),
            (ReactionRole::Agent, &mut reaction.agents),
            (ReactionRole::Product, &mut reaction.products),
        ] {
            for (template, graph) in templates.iter_mut().enumerate() {
                let parsed = cosmolkit_cx::parse_cx_extensions_with_atom_window(
                    cx,
                    start_atom as usize,
                    graph.num_atoms() as u32 as usize,
                )
                .map_err(|source| ReactionParseError::CxParse {
                    role,
                    template,
                    start_atom: start_atom as usize,
                    start_bond: start_bond as usize,
                    source,
                })?;
                cosmolkit_search::apply_cx_to_query_graph_with_offsets(
                    graph,
                    &parsed,
                    start_atom as usize,
                    start_bond as usize,
                )
                .map_err(|source| ReactionParseError::CxLowering {
                    role,
                    template,
                    start_atom: start_atom as usize,
                    start_bond: start_bond as usize,
                    source,
                })?;
                start_atom = start_atom.wrapping_add(graph.num_atoms() as u32);
                start_bond = start_bond.wrapping_add(graph.num_bonds() as u32);
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
    reaction.set_implicit_properties_source(true);
    Ok(reaction)
}

fn reaction_from_smiles_source(
    text: &str,
    params: &ReactionParseParams,
) -> Result<Reaction, ReactionParseError> {
    // RDKit❗✔️: std::unique_ptr<ChemicalReaction> ReactionFromSmiles(
    // RDKit❗✔️:     const std::string &origText, const ReactionSmartsParserParams &options) {
    // RDKit❗✔️:   return parseReaction(origText, options, true);
    // RDKit❗✔️: }
    parse_reaction_source(text, params, true)
}

fn reaction_from_smarts_source(
    text: &str,
    params: &ReactionParseParams,
) -> Result<Reaction, ReactionParseError> {
    // RDKit❗✔️: std::unique_ptr<ChemicalReaction> ReactionFromSmarts(
    // RDKit❗✔️:     const std::string &origText, const ReactionSmartsParserParams &options) {
    // RDKit❗✔️:   return parseReaction(origText, options, false);
    // RDKit❗✔️: }
    parse_reaction_source(text, params, false)
}

pub fn parse_smirks(text: &str) -> Result<Reaction, ReactionParseError> {
    reaction_from_smarts_source(text, &ReactionParseParams::default())
}

#[cfg(test)]
mod complete_split_components_source_tests {
    use super::*;
    #[test]
    fn grouped_dots_empty_components_and_source_trailing_rules() {
        for (input, expected) in [
            ("", vec![]),
            (".", vec![""]),
            ("..", vec!["", ""]),
            (".A.", vec!["", "A"]),
            ("A..B", vec!["A", "", "B"]),
            ("(A.B).C", vec!["A.B", "C"]),
            ("A(B.C).D", vec!["A(B.C)", "D"]),
            ("((A.B)).C", vec!["(A.B)", "C"]),
            ("(A)B.C", vec!["A)", "C"]),
            ("()..", vec!["", ""]),
        ] {
            let expected: Vec<&[u8]> = expected.iter().map(|s| s.as_bytes()).collect();
            assert_eq!(
                split_components_source(input.as_bytes()).unwrap(),
                expected,
                "{input:?}"
            );
        }
    }
    #[test]
    fn unsigned_parenthesis_underflow_and_closed_block_state_are_not_repaired() {
        assert_eq!(
            split_components_source(b")(A.B").unwrap(),
            vec![b")(A".as_slice(), b"B".as_slice()]
        );
        assert_eq!(
            split_components_source(b"(A))B.C").unwrap(),
            vec![b"A))B.".as_slice()]
        );
        assert_eq!(
            split_components_source(b"(").unwrap(),
            vec![b"(".as_slice()]
        );
        assert_eq!(
            split_components_source(b")(").unwrap(),
            vec![b")(".as_slice()]
        );
    }
    #[test]
    fn counted_nul_and_non_utf8_component_bytes_survive_exact_source_substrings() {
        let input = b"(A\0\xff).B";
        let parts = split_components_source(input).unwrap();
        assert_eq!(parts, vec![b"A\0\xff".as_slice(), b"B".as_slice()]);
        assert_eq!(parts[0].as_ptr(), input[1..].as_ptr());
        assert_eq!(
            split_components_source(b"\xff.\0").unwrap(),
            vec![b"\xff".as_slice(), b"\0".as_slice()]
        );
        let unicode = "(A)é.C";
        assert_eq!(
            split_components_source(unicode.as_bytes()).unwrap(),
            vec![b"A)\xc3".as_slice(), b"C".as_slice()]
        );
    }
    #[test]
    fn native_substr_count_clamp_and_existing_checked_utf8_projection_are_explicit() {
        assert_eq!(component_substr_source(b"abc", 1, u32::MAX).unwrap(), b"bc");
        assert_eq!(component_substr_source(b"abc", 3, u32::MAX).unwrap(), b"");
        assert!(matches!(
            component_substr_source(b"abc", 4, 0),
            Err(ReactionParseError::ComponentBounds { start: 4, end: 4 })
        ));
        assert_eq!(split_components("(A.B).C").unwrap(), vec!["A.B", "C"]);
        assert!(matches!(
            split_components("(A)é.C"),
            Err(ReactionParseError::ComponentBounds { start: 1, end: 4 })
        ));
    }
}

#[cfg(test)]
mod complete_construct_component_source_tests {
    use super::*;
    use cosmolkit_model::{
        Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, PropertyValue,
    };
    fn params(use_smiles: bool) -> ReactionParseParams {
        ReactionParseParams {
            use_smiles,
            ..Default::default()
        }
    }
    #[test]
    fn dispatch_replacements_and_native_query_presence_follow_selected_parser() {
        for use_smiles in [false, true] {
            let mut params = params(use_smiles);
            params.replacements.insert("{Q}".into(), "[C:7]O".into());
            let graph = construct_component(
                "{Q}",
                &params,
                (&params).use_smiles,
                ReactionRole::Reactant,
                3,
            )
            .unwrap();
            assert_eq!(graph.num_atoms(), 2);
            assert_eq!(graph.num_bonds(), 1);
            assert_eq!(graph.atoms()[0].atom_map(), Some(7));
            assert!(
                graph
                    .atoms()
                    .iter()
                    .all(|atom| atom.predicate_is_carrier_derived() == use_smiles)
            );
            assert!(
                graph
                    .bonds()
                    .iter()
                    .all(|bond| bond.predicate_is_carrier_derived() == use_smiles)
            );
        }
    }
    #[test]
    fn source_disables_component_name_and_cx_parsing_without_fallback_or_context_loss() {
        for use_smiles in [false, true] {
            for text in ["C trailing_name", "C |$tag$|"] {
                let error = construct_component(
                    text,
                    &params(use_smiles),
                    (&params(use_smiles)).use_smiles,
                    ReactionRole::Product,
                    6,
                )
                .unwrap_err();
                match error {
                    ReactionParseError::Smarts {
                        role: ReactionRole::Product,
                        template: 6,
                        text: actual,
                        ..
                    } if !use_smiles => assert_eq!(actual.as_bytes(), text.as_bytes()),
                    ReactionParseError::Smiles {
                        role: ReactionRole::Product,
                        template: 6,
                        text: actual,
                        ..
                    } if use_smiles => assert_eq!(actual.as_bytes(), text.as_bytes()),
                    other => panic!("wrong parser/context: {other:?}"),
                }
            }
        }
    }
    #[test]
    fn source_keeps_explicit_hydrogens_in_both_modes() {
        for use_smiles in [false, true] {
            let graph = construct_component(
                "[H][C:7]",
                &params(use_smiles),
                (&params(use_smiles)).use_smiles,
                ReactionRole::Agent,
                0,
            )
            .unwrap();
            assert_eq!(graph.num_atoms(), 2);
            assert_eq!(graph.num_bonds(), 1);
            assert_eq!(graph.atoms()[0].atomic_number(), 1);
            assert_eq!(graph.atoms()[1].atom_map(), Some(7));
        }
    }
    #[test]
    fn sanitize_flag_is_forwarded_only_to_smiles() {
        let text = "C(C)(C)(C)(C)C";
        let mut smarts = params(false);
        smarts.sanitize = true;
        assert_eq!(
            construct_component(
                text,
                &smarts,
                (&smarts).use_smiles,
                ReactionRole::Reactant,
                0
            )
            .unwrap()
            .num_atoms(),
            6
        );
        let smiles = params(true);
        assert_eq!(
            construct_component(
                text,
                &smiles,
                (&smiles).use_smiles,
                ReactionRole::Reactant,
                0
            )
            .unwrap()
            .num_atoms(),
            6
        );
        let mut sanitize = smiles;
        sanitize.sanitize = true;
        assert!(matches!(
            construct_component(
                text,
                &sanitize,
                (&sanitize).use_smiles,
                ReactionRole::Reactant,
                0
            ),
            Err(ReactionParseError::Smiles { .. })
        ));
    }
    #[test]
    fn detached_carrier_projection_preserves_ordered_byte_properties_and_mixed_conformer_order() {
        let mut record = cosmolkit_smiles::parse_smiles_complete_source(
            "[C:7]O",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                skip_cleanup: true,
                allow_cxsmiles: false,
                parse_name: false,
                ..Default::default()
            },
        )
        .unwrap();
        record
            .properties
            .set_prop("z", PropertyValue::Int(7))
            .unwrap();
        record
            .properties
            .set_prop("a", PropertyValue::Bool(false))
            .unwrap();
        record
            .properties
            .set_prop(vec![0, 255], PropertyValue::UInt(0))
            .unwrap();
        let properties: Vec<_> = record
            .properties
            .ordered_props()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        record.coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(11, vec![[1.0, 2.0], [3.0, 4.0]])],
            conformers_3d: vec![Conformer3D::new(
                4,
                vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                true,
            )],
            source_conformer_order: Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
            ..Default::default()
        };
        let graph = project_smiles_component_record(record, ReactionRole::Product, 8).unwrap();
        assert_eq!(
            graph
                .ordered_props()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect::<Vec<_>>(),
            properties
        );
        assert_eq!(
            graph.source_conformer_order(),
            Some([CoordinateDimension::ThreeD, CoordinateDimension::TwoD].as_slice())
        );
        let coordinates = graph.coordinate_block(None);
        assert_eq!(coordinates.conformers_2d[0].id(), 11);
        assert_eq!(coordinates.conformers_3d[0].id(), 4);
        assert!(
            graph
                .atoms()
                .iter()
                .all(|atom| atom.predicate_is_carrier_derived())
        );
        assert!(
            graph
                .bonds()
                .iter()
                .all(|bond| bond.predicate_is_carrier_derived())
        );
    }
}

#[cfg(test)]
mod complete_remove_spaces_source_tests {
    use super::*;

    #[test]
    fn removes_only_adjacent_space_and_tab_in_source_order() {
        for (input, pos, expected) in [
            ("A \t>\t B", 3, "A>B"),
            (" \t> \t", 2, ">"),
            ("> \tA", 0, ">A"),
            ("A \t>", 3, "A>"),
            ("A\n > \rB", 3, "A\n>\rB"),
            ("é \t>\t λ", 4, "é>λ"),
            ("A\0 > \tB", 3, "A\0>B"),
        ] {
            let mut text = input.to_owned();
            let capacity = text.capacity();
            remove_spaces_around(&mut text, pos).unwrap();
            assert_eq!(text, expected);
            assert_eq!(text.capacity(), capacity);
        }
    }

    #[test]
    fn source_endpoint_nul_is_defined_even_one_past_delimiter_position() {
        for input in ["", "A", "A \t", "é"] {
            let mut text = input.to_owned();
            remove_spaces_around(&mut text, input.len() + 1).unwrap();
            assert_eq!(text, input);
        }
        let mut text = "A \t".to_owned();
        remove_spaces_around(&mut text, 3).unwrap();
        assert_eq!(text, "A");
    }

    #[test]
    fn undefined_source_index_is_structured_and_keeps_prior_erase_effects() {
        let mut text = " \tA".to_owned();
        assert!(matches!(remove_spaces_around(&mut text, usize::MAX),
            Err(ReactionParseError::ComponentBounds { start, end })
                if start == usize::MAX - 1 && end == start));
        assert_eq!(text, "A");
        let mut text = "é".to_owned();
        assert!(matches!(
            remove_spaces_around(&mut text, 4),
            Err(ReactionParseError::ComponentBounds { start: 3, end: 3 })
        ));
        assert_eq!(text, "é");
    }
}

#[cfg(test)]
mod complete_parse_reaction_source_tests {
    use super::*;
    use cosmolkit_model::PropertyValue;

    #[test]
    fn c_locale_edge_trim_includes_vertical_tab_and_preserves_unicode_space() {
        assert_eq!(trim_source_text("\x0b\r\n\t C>>N \x0c\x0b"), "C>>N");
        assert_eq!(trim_source_text("\u{a0}C>>N\u{a0}"), "\u{a0}C>>N\u{a0}");
        let reaction = parse_smirks("\x0b C>>N \x0b").unwrap();
        assert_eq!(reaction.num_reactant_templates(), 1);
        assert_eq!(reaction.num_product_templates(), 1);
    }

    #[test]
    fn name_truncation_refinds_separators_and_dative_arrows_are_not_separators() {
        let reaction = parse_smirks("C>>N reaction > name").unwrap();
        assert_eq!(reaction.num_product_templates(), 1);
        assert!(matches!(
            parse_smirks("C name >>N"),
            Err(ReactionParseError::Separators { count: 0 })
        ));
        assert!(matches!(
            parse_smirks("C>>N>>O"),
            Err(ReactionParseError::MultiStep { count: 4 })
        ));
        assert!(matches!(
            parse_smirks("C->N"),
            Err(ReactionParseError::Separators { count: 0 })
        ));
        let reaction = parse_smirks("C->N>>C->N").unwrap();
        assert_eq!(reaction.reactants[0].num_atoms(), 2);
        assert_eq!(reaction.products[0].num_bonds(), 1);
    }

    #[test]
    fn grouping_dots_and_adjacent_whitespace_keep_source_component_order() {
        let reaction = parse_smirks(" (C . O) >> (N . C) ").unwrap();
        assert_eq!(reaction.num_reactant_templates(), 1);
        assert_eq!(reaction.reactants[0].num_atoms(), 2);
        assert_eq!(reaction.products[0].num_atoms(), 2);
        let reaction = parse_smirks(".C..>>N.").unwrap();
        assert_eq!(
            reaction
                .reactants
                .iter()
                .map(QueryGraph::num_atoms)
                .collect::<Vec<_>>(),
            [0, 1, 0]
        );
        assert_eq!(
            reaction
                .products
                .iter()
                .map(QueryGraph::num_atoms)
                .collect::<Vec<_>>(),
            [1]
        );
    }

    #[test]
    fn component_byte_slices_reach_smarts_owner_without_utf8_repair() {
        match parse_smirks("(C)é>>N").unwrap_err() {
            ReactionParseError::Smarts {
                role: ReactionRole::Reactant,
                template: 0,
                text,
                ..
            } => assert_eq!(text.as_bytes(), &[b'C', b')', 0xc3]),
            other => panic!("source byte component was not forwarded: {other:?}"),
        }
        let error = construct_component(
            [b'C', 0xff],
            &ReactionParseParams {
                use_smiles: true,
                ..Default::default()
            },
            (&ReactionParseParams {
                use_smiles: true,
                ..Default::default()
            })
                .use_smiles,
            ReactionRole::Agent,
            2,
        )
        .unwrap_err();
        assert!(
            matches!(error, ReactionParseError::ComponentEncoding { role: ReactionRole::Agent, template: 2, text, .. }
            if text.as_bytes() == [b'C', 0xff])
        );
    }

    #[test]
    fn source_role_failure_order_is_reactants_then_products_then_agents() {
        assert!(matches!(
            parse_smirks("[>[>["),
            Err(ReactionParseError::Smarts {
                role: ReactionRole::Reactant,
                template: 0,
                ..
            })
        ));
        assert!(matches!(
            parse_smirks("C>[>["),
            Err(ReactionParseError::Smarts {
                role: ReactionRole::Product,
                template: 0,
                ..
            })
        ));
        assert!(matches!(
            parse_smirks("C>[>C"),
            Err(ReactionParseError::Smarts {
                role: ReactionRole::Agent,
                template: 0,
                ..
            })
        ));
    }

    #[test]
    fn agents_fragment_before_global_cx_reactant_agent_product_windows() {
        let reaction = parse_smirks("C>O.N>C |$react;agent1;agent2;product$|").unwrap();
        assert_eq!(reaction.num_agent_templates(), 2);
        for (graph, label) in [
            (&reaction.reactants[0], "react"),
            (&reaction.agents[0], "agent1"),
            (&reaction.agents[1], "agent2"),
            (&reaction.products[0], "product"),
        ] {
            assert_eq!(
                graph.atoms()[0].prop("atomLabel"),
                Some(&PropertyValue::from(label))
            );
            assert!(graph.prop("_cxsmilesLabelsProcessed").is_none());
            assert!(graph.prop("_cxsmiles_sgroup_tracker").is_none());
        }
    }

    #[test]
    fn source_cx_guard_ignores_strict_flag_and_empty_template_loops_skip_invalid_cx() {
        for strict in [false, true] {
            let params = ReactionParseParams {
                strict_cxsmiles: strict,
                ..Default::default()
            };
            assert!(matches!(
                parse_smirks_with_params("C>>C |invalid|", &params),
                Err(ReactionParseError::CxParse {
                    role: ReactionRole::Reactant,
                    ..
                })
            ));
        }
        let disabled = ReactionParseParams {
            allow_cxsmiles: false,
            ..Default::default()
        };
        assert!(parse_smirks_with_params("C>>C |invalid|", &disabled).is_ok());
        let reaction = parse_smirks(">> |invalid|").unwrap();
        assert_eq!(reaction.num_reactant_templates(), 0);
        assert!(reaction.implicit_properties);
    }

    #[test]
    fn both_parser_modes_finish_cleanup_then_enable_implicit_properties_without_init() {
        for use_smiles in [false, true] {
            let params = ReactionParseParams {
                use_smiles,
                ..Default::default()
            };
            let reaction = parse_smirks_with_params("[H][C:1]>O.[H]>[H][C:1]", &params).unwrap();
            assert!(reaction.implicit_properties);
            assert!(reaction.needs_init);
            assert_eq!(reaction.reactants[0].num_atoms(), 2);
            assert_eq!(reaction.num_agent_templates(), 2);
            for graph in reaction
                .reactants
                .iter()
                .chain(&reaction.agents)
                .chain(&reaction.products)
            {
                assert!(
                    graph
                        .atoms()
                        .iter()
                        .all(|atom| atom.prop("_SmilesStart").is_none())
                );
                assert!(
                    graph
                        .bonds()
                        .iter()
                        .all(|bond| bond.bond().prop("_cxsmilesBondIdx").is_none())
                );
            }
        }
    }
}

#[cfg(test)]
mod complete_reaction_from_smarts_source_tests {
    use super::*;

    #[test]
    fn source_smarts_selector_is_false_independent_of_projection_mode_field() {
        let mut params = ReactionParseParams {
            use_smiles: true,
            sanitize: true,
            ..Default::default()
        };
        params
            .replacements
            .insert("{Q}".into(), "C(C)(C)(C)(C)C".into());
        let reaction = reaction_from_smarts_source("{Q}>>C", &params).unwrap();
        assert_eq!(reaction.reactants[0].num_atoms(), 6);
        assert!(
            reaction.reactants[0]
                .atoms()
                .iter()
                .all(|a| !a.predicate_is_carrier_derived())
        );
        assert!(matches!(
            parse_smirks_with_params("{Q}>>C", &params),
            Err(ReactionParseError::Smiles {
                role: ReactionRole::Reactant,
                ..
            })
        ));
    }

    #[test]
    fn full_options_reference_is_forwarded_to_the_one_parser_body() {
        let mut params = ReactionParseParams {
            allow_cxsmiles: false,
            ..Default::default()
        };
        params.replacements.insert("{Q}".into(), "[N:7]".into());
        let reaction = reaction_from_smarts_source("{Q}>>{Q} |invalid|", &params).unwrap();
        assert_eq!(reaction.reactants[0].atoms()[0].atom_map(), Some(7));
        assert_eq!(reaction.products[0].atoms()[0].mol_inversion_flag(), None);
        assert!(reaction.implicit_properties);
        assert!(reaction.needs_init);
        params.allow_cxsmiles = true;
        assert!(matches!(
            reaction_from_smarts_source("{Q}>>{Q} |invalid|", &params),
            Err(ReactionParseError::CxParse {
                role: ReactionRole::Reactant,
                ..
            })
        ));
    }
}

#[cfg(test)]
mod complete_reaction_from_smiles_source_tests {
    use super::*;

    #[test]
    fn native_smiles_selector_is_true_for_both_combined_projection_mode_values() {
        for mode in [false, true] {
            let params = ReactionParseParams {
                use_smiles: mode,
                ..Default::default()
            };
            let reaction = reaction_from_smiles_source("[H][C:7]>>[H][C:7]", &params).unwrap();
            assert_eq!(reaction.reactants[0].num_atoms(), 2);
            assert!(
                reaction.reactants[0]
                    .atoms()
                    .iter()
                    .all(|a| a.predicate_is_carrier_derived())
            );
            assert!(
                reaction.products[0]
                    .bonds()
                    .iter()
                    .all(|b| b.predicate_is_carrier_derived())
            );
            assert!(reaction.implicit_properties);
            assert!(reaction.needs_init);
        }
    }

    #[test]
    fn native_smiles_overload_forwards_replacements_and_sanitize_options() {
        let mut params = ReactionParseParams::default();
        params
            .replacements
            .insert("{Q}".into(), "C(C)(C)(C)(C)C".into());
        assert_eq!(
            reaction_from_smiles_source("{Q}>>C", &params)
                .unwrap()
                .reactants[0]
                .num_atoms(),
            6
        );
        params.sanitize = true;
        assert!(matches!(
            reaction_from_smiles_source("{Q}>>C", &params),
            Err(ReactionParseError::Smiles {
                role: ReactionRole::Reactant,
                template: 0,
                ..
            })
        ));
        assert!(reaction_from_smarts_source("{Q}>>C", &params).is_ok());
    }
}
