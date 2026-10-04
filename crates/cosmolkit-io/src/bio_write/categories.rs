//! Block-name and `_entry.id` preamble for the coordinate writer.
//!
//! This module ports only the pinned preamble of Gemmi's
//! `update_mmcif_block` (to_mmcif.cpp:460-473) under this packet's fixed
//! profile: `block_name=true` and `entry=true`, every other group off. No
//! existing-document update API exists here and no additional group switch
//! is modeled.

use cosmolkit_bio::BioStructureData;

use crate::cif::CifBlock;
use crate::cif::quote_cif_value;

/// Gemmi `is_valid_block_name` (to_mmcif.cpp:305-307).
fn is_valid_block_name(name: &str) -> bool {
    // Gemmi✔️✔️: bool is_valid_block_name(const std::string& name) {
    // Gemmi✔️✔️:   return !name.empty() &&
    // Gemmi✔️✔️:          std::all_of(name.begin(), name.end(), [](char c){ return c >= '!' && c <= '~'; });
    // Gemmi✔️✔️: }
    // Behavior: nonempty and every byte printable ASCII ('!'..'~'); the
    // CK name is a Rust String, so a non-ASCII name simply fails the range
    // test exactly like the source's char comparison.
    // Complexity: single linear byte scan, no allocation.
    !name.is_empty() && name.bytes().all(|b| (b'!'..=b'~').contains(&b))
}

/// Gemmi `update_mmcif_block` preamble for a fresh block under the fixed
/// enabled `block_name`/`entry` profile; the source zero-model early return
/// is preserved and later group sections are not modeled in this scope.
pub(crate) fn update_mmcif_block_preamble(data: &BioStructureData, block: &mut CifBlock) {
    // Gemmi✔️✔️: void update_mmcif_block(const Structure& st, cif::Block& block, MmcifOutputGroups groups) {
    // Gemmi✔️✔️:   if (st.models.empty())
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:   if (groups.block_name)
    // Gemmi✔️✔️:     block.name = is_valid_block_name(st.name) ? st.name : "model";
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:   auto e_id = st.info.find("_entry.id");
    // Gemmi✔️✔️:   std::string id = cif::quote(e_id != st.info.end() ? e_id->second : block.name);
    // Gemmi✔️✔️:   if (groups.entry)
    // Gemmi✔️✔️:     block.set_pair("_entry.id", id);
    // Gemmi✔️✔️:   else if (const std::string* val = block.find_value("_entry.id"))
    // Gemmi✔️✔️:     id = *val;
    // Behavior: with no models the block is returned untouched (the same
    // early return guards the whole source function). The block name is the
    // structure name when it is a valid CIF block name, otherwise the literal
    // "model"; the entry id prefers the retained `_entry.id` info value and
    // falls back to the just-assigned block name, and the quoted id replaces
    // or inserts the `_entry.id` pair through the whole-block span, exactly
    // like `block.set_pair`. The `else` fallback branch is unreachable under
    // this packet's permanently enabled `entry` profile and is not modeled;
    // `database_status` and every later group section are out of scope.
    // Complexity: one models emptiness check, one name scan, one BTreeMap
    // lookup and the set_pair linear scan with at most one insertion.
    if data.models().is_empty() {
        return;
    }

    let state = data.source_state();
    block.set_name(if is_valid_block_name(&state.name) {
        state.name.clone()
    } else {
        "model".to_string()
    });

    let id = match state.info.get("_entry.id") {
        Some(value) => quote_cif_value(value.clone()),
        None => quote_cif_value(block.name().to_string()),
    };
    block.set_pair_in_category(None, "_entry.id", id);
}

#[cfg(test)]
mod tests {
    use super::update_mmcif_block_preamble;
    use crate::cif::{CifBlock, CifCheckLevel, read_cif_document};
    use cosmolkit_bio::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureParts,
        BioStructureSourceState, ChainKind, ChainSourceIds, EntityKind, ResidueInfoKind,
        ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;

    fn single_model_parts(name: &str, entry_id: Option<&str>) -> BioStructureParts {
        let mut parts = BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            ..empty_base()
        };
        let mut info = BTreeMap::new();
        if let Some(id) = entry_id {
            info.insert("_entry.id".to_string(), id.to_string());
        }
        parts.source_state = BioStructureSourceState {
            name: name.to_string(),
            info,
            ..BioStructureSourceState::default()
        };
        parts.models = vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))];
        parts.chains = vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::default(),
        )];
        parts.residues = vec![BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        )];
        parts.atoms = vec![BioAtomRow::new(
            BioResidueId::new(0),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            AtomSourceIds::default(),
        )];
        parts.coordinates = BioCoordinateBlock::new(vec![[0.0, 0.0, 0.0]]);
        parts
    }

    fn empty_base() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: Vec::new(),
            chains: Vec::new(),
            residues: Vec::new(),
            atoms: Vec::new(),
            entities: Vec::new(),
            connections: Vec::new(),
            cispeps: Vec::new(),
            mod_residues: Vec::new(),
            helices: Vec::new(),
            sheets: Vec::new(),
            metadata: cosmolkit_bio::BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::new(Vec::new()),
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        }
    }

    fn fresh_block() -> CifBlock {
        read_cif_document("data_fresh\n", "p01", CifCheckLevel::Syntax)
            .unwrap()
            .blocks()[0]
            .clone()
    }

    #[test]
    fn bio_pdbscope_p01_valid_name_and_entry_selection_from_metadata() {
        // to_mmcif.cpp:460-473: entry id prefers st.info["_entry.id"], falls
        // back to the just-assigned block name; both spellings below derive
        // from cif::quote's rules, not CK output.
        let with_entry =
            BioStructureData::from_parts(single_model_parts("1abc", Some("1ABC"))).unwrap();
        let mut block = fresh_block();
        update_mmcif_block_preamble(&with_entry, &mut block);
        assert_eq!(block.name(), "1abc");
        assert_eq!(block.find_value("_entry.id").unwrap().raw(), "1ABC");

        let without_entry = BioStructureData::from_parts(single_model_parts("1abc", None)).unwrap();
        let mut block = fresh_block();
        update_mmcif_block_preamble(&without_entry, &mut block);
        assert_eq!(block.name(), "1abc");
        assert_eq!(block.find_value("_entry.id").unwrap().raw(), "1abc");

        // Quoting rules: a space forces quotes; an empty info value keeps
        // the opening/closing quote pair (''), never a dot or question mark.
        let spaced = BioStructureData::from_parts(single_model_parts("n", Some("1A B"))).unwrap();
        let mut block = fresh_block();
        update_mmcif_block_preamble(&spaced, &mut block);
        assert_eq!(block.find_value("_entry.id").unwrap().raw(), "'1A B'");
        let empty_id = BioStructureData::from_parts(single_model_parts("n", Some(""))).unwrap();
        let mut block = fresh_block();
        update_mmcif_block_preamble(&empty_id, &mut block);
        assert_eq!(block.find_value("_entry.id").unwrap().raw(), "''");
    }

    #[test]
    fn bio_pdbscope_p01_invalid_and_empty_names_fall_back_to_model() {
        // is_valid_block_name (to_mmcif.cpp:305-307): nonempty and every
        // byte in '!'..'~'. A space, a newline and a non-ASCII name each
        // fail; the fallback is the literal "model".
        for invalid in ["", "a b", "a\nb", "a\u{e9}b"] {
            let data = BioStructureData::from_parts(single_model_parts(invalid, None)).unwrap();
            let mut block = fresh_block();
            update_mmcif_block_preamble(&data, &mut block);
            assert_eq!(block.name(), "model", "name {invalid:?} must fall back");
            assert_eq!(block.find_value("_entry.id").unwrap().raw(), "model");
        }
    }

    #[test]
    fn bio_pdbscope_p01_no_models_leaves_block_untouched() {
        // The source early return guards the entire function: a structure
        // without models writes neither the name nor the entry pair.
        let empty = BioStructureData::from_parts(empty_base()).unwrap();
        let mut block = fresh_block();
        update_mmcif_block_preamble(&empty, &mut block);
        assert_eq!(block.name(), "fresh", "parsed block name is preserved");
        assert!(block.find_value("_entry.id").is_none());
        assert!(block.items().is_empty());
    }
}
