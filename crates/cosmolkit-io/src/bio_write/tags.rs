//! Source enum/tag conversions and the atom-row entity-id lookup.

use cosmolkit_bio::{BioCalcFlag, BioEntityRow, BioResidueRow, EntityKind, find_residue_info};

use super::BioMmcifWriteError;
use crate::cif::quote_cif_value;
/// Gemmi `use_hetatm`: decide the `_atom_site.group_PDB` record tag.
pub(crate) fn use_hetatm(residue: &BioResidueRow) -> bool {
    // Gemmi✔️✔️: bool use_hetatm(const Residue& res) {
    // Gemmi✔️✔️:   if (res.het_flag == 'H')
    // Gemmi✔️✔️:     return true;
    // Gemmi✔️✔️:   if (res.het_flag == 'A')
    // Gemmi✔️✔️:     return false;
    // Gemmi✔️✔️:   if (res.entity_type == EntityType::Branched ||
    // Gemmi✔️✔️:       res.entity_type == EntityType::NonPolymer ||
    // Gemmi✔️✔️:       res.entity_type == EntityType::Water)
    // Gemmi✔️✔️:     return true;
    // Gemmi✔️✔️:   return !find_tabulated_residue(res.name).is_standard();
    // Gemmi✔️✔️: }
    // Behavior: the source het_flag byte drives the first two branches
    // (BIO stores it as Option<u8>; None falls through exactly like any
    // byte other than 'H'/'A'), then the entity kind, then the tabulated
    // residue table decides. No kind is guessed or reclassified.
    // Complexity: O(1) except the final table lookup, which is the same
    // tabulated lookup the source performs.
    match residue.het_flag() {
        Some(b'H') => return true,
        Some(b'A') => return false,
        _ => {}
    }
    matches!(
        residue.entity_kind(),
        EntityKind::Branched | EntityKind::NonPolymer | EntityKind::Water
    ) || !find_residue_info(residue.name().as_str()).is_standard()
}

/// Gemmi `_atom_site.calc_flag` text: the `".\0.\0d\0c\0dum"` lookup table.
pub(crate) const fn calc_flag_text(flag: BioCalcFlag) -> &'static str {
    // Gemmi✔️✔️: if (has_calc_flag)
    // Gemmi✔️✔️:   vv.emplace_back(&".\0.\0d\0c\0dum"[2 * (int) atom.calc_flag]);
    // Behavior: the source indexes a NUL-separated literal by twice the
    // enum ordinal: NotSet and NoHydrogen both print ".", Determined "d",
    // Calculated "c", Dummy "dum". The variant order is part of the source
    // ABI and is mirrored by BioCalcFlag's declaration order.
    // Complexity: O(1) table selection.
    match flag {
        BioCalcFlag::NotSet => ".",
        BioCalcFlag::NoHydrogen => ".",
        BioCalcFlag::Determined => "d",
        BioCalcFlag::Calculated => "c",
        BioCalcFlag::Dummy => "dum",
    }
}

/// Gemmi atom-row entity id: the owning entity's quoted name, or the
/// residue's own label entity id as `string_or_dot`.
pub(crate) fn entity_id_for_residue(
    residue: &BioResidueRow,
    entities: &[BioEntityRow],
) -> Result<String, BioMmcifWriteError> {
    // Gemmi✔️✔️: inline Entity* find_entity_of_subchain(const std::string& subchain_id,
    // Gemmi✔️✔️:                                        std::vector<Entity>& entities) {
    // Gemmi✔️✔️:   if (!subchain_id.empty())
    // Gemmi✔️✔️:     for (Entity& ent : entities)
    // Gemmi✔️✔️:       if (in_vector(subchain_id, ent.subchains))
    // Gemmi✔️✔️:         return &ent;
    // Gemmi✔️✔️:   return nullptr;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: if (const Entity* ent = gemmi::find_entity_of_subchain(res.subchain, st.entities))
    // Gemmi✔️✔️:   entity_id = cif::quote(ent->name);
    // Gemmi✔️✔️: else
    // Gemmi✔️✔️:   entity_id = string_or_dot(res.entity_id);
    // Behavior: an empty/absent subchain never matches; the first entity
    // listing the subchain wins; its source entity name is quoted. Without
    // a matching entity the residue's own label entity id follows the
    // string_or_dot rule.
    // Complexity: O(entities x subchains) linear scan like the source, plus
    // one quoting per call; no map or index structure is added.
    let subchain = residue.source().subchain_id().unwrap_or_default();
    if !subchain.is_empty() {
        for entity in entities {
            if entity
                .subchains()
                .iter()
                .any(|candidate| candidate == subchain)
            {
                return Ok(quote_cif_value(
                    entity.source().source_entity_id().to_string(),
                ));
            }
        }
    }
    let fallback = residue.source().label_entity_id().map_or_else(
        || ".".to_string(),
        |id| {
            if id.is_empty() {
                ".".to_string()
            } else {
                quote_cif_value(id.to_string())
            }
        },
    );
    Ok(fallback)
}

#[cfg(test)]
mod tests {
    use super::{BioMmcifWriteError, calc_flag_text, entity_id_for_residue, use_hetatm};
    use cosmolkit_bio::{
        BioCalcFlag, BioChainId, BioEntityRow, BioResidueRow, BioRowSpan, BioSiftsUnpResidue,
        BioStructureParts, ChainKind, ChainSourceIds, EntityKind, EntitySourceIds, PolymerKind,
        ResidueInfoKind, ResidueName, ResidueSourceIds, find_residue_info,
    };

    fn residue(
        name: &str,
        entity_kind: EntityKind,
        het_flag: Option<u8>,
        subchain: Option<&str>,
        label_entity_id: Option<&str>,
    ) -> BioResidueRow {
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 0).unwrap(),
            ResidueName::from_ascii(name.as_bytes()).unwrap(),
            ResidueInfoKind::Aa,
            entity_kind,
            None,
            het_flag,
            ResidueSourceIds::new(
                None,
                None,
                None,
                subchain.map(str::to_string),
                label_entity_id.map(str::to_string),
            )
            .unwrap(),
            BioSiftsUnpResidue::default(),
        )
    }

    fn entity(name: &str, kind: EntityKind, subchains: &[&str]) -> BioEntityRow {
        BioEntityRow::new(
            kind,
            PolymerKind::Unknown,
            false,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            subchains.iter().map(|value| value.to_string()).collect(),
            EntitySourceIds::new(name.to_string()),
        )
    }

    #[test]
    fn bio_pdbscope_v04_het_flag_entity_kind_and_tabulated_fallback() {
        // Explicit flags short-circuit everything.
        assert!(use_hetatm(&residue(
            "ALA",
            EntityKind::Unknown,
            Some(b'H'),
            None,
            None
        )));
        assert!(!use_hetatm(&residue(
            "HOH",
            EntityKind::Water,
            Some(b'A'),
            None,
            None
        )));
        // Entity kinds force HETATM before the table is consulted.
        for kind in [
            EntityKind::Branched,
            EntityKind::NonPolymer,
            EntityKind::Water,
        ] {
            assert!(
                use_hetatm(&residue("ALA", kind, None, None, None)),
                "{kind:?}"
            );
        }
        // Fallback is exactly the tabulated standard-ness of the name.
        assert!(!use_hetatm(&residue(
            "ALA",
            EntityKind::Polymer,
            None,
            None,
            None
        )));
        assert!(use_hetatm(&residue(
            "HOH",
            EntityKind::Polymer,
            None,
            None,
            None
        )));
        // Unknown het_flag bytes fall through like the source comparison.
        assert!(!use_hetatm(&residue(
            "ALA",
            EntityKind::Polymer,
            Some(b'?'),
            None,
            None
        )));
        // Any other byte than 'H'/'A' does not short-circuit.
        assert!(use_hetatm(&residue(
            "HOH",
            EntityKind::Polymer,
            Some(b'X'),
            None,
            None
        )));
        // The tabulated lookup used above is the shared owner.
        assert!(find_residue_info("ALA").is_standard());
        assert!(!find_residue_info("HOH").is_standard());
    }

    #[test]
    fn bio_pdbscope_v04_calc_flag_table_is_source_exact() {
        assert_eq!(calc_flag_text(BioCalcFlag::NotSet), ".");
        assert_eq!(calc_flag_text(BioCalcFlag::NoHydrogen), ".");
        assert_eq!(calc_flag_text(BioCalcFlag::Determined), "d");
        assert_eq!(calc_flag_text(BioCalcFlag::Calculated), "c");
        assert_eq!(calc_flag_text(BioCalcFlag::Dummy), "dum");
    }

    #[test]
    fn bio_pdbscope_v04_entity_lookup_order_and_fallback() -> Result<(), BioMmcifWriteError> {
        let polymer = entity("polymerEnt", EntityKind::Polymer, &["A", "B"]);
        let water = entity("waterEnt", EntityKind::Water, &["C"]);
        let entities = vec![polymer, water];

        // Subchain hit quotes the first matching entity's name.
        let in_a = residue("ALA", EntityKind::Polymer, None, Some("A"), None);
        assert_eq!(entity_id_for_residue(&in_a, &entities)?, "polymerEnt");
        let in_c = residue("HOH", EntityKind::Water, None, Some("C"), None);
        assert_eq!(entity_id_for_residue(&in_c, &entities)?, "waterEnt");

        // Empty or absent subchain never matches, even if an entity lists "".
        let empty_listed = entity("emptyEnt", EntityKind::NonPolymer, &[""]);
        let no_subchain = residue("ALA", EntityKind::Polymer, None, None, Some("labelEnt"));
        assert_eq!(
            entity_id_for_residue(&no_subchain, &[empty_listed])?,
            "labelEnt"
        );
        let empty_subchain = residue("ALA", EntityKind::Polymer, None, Some(""), None);
        assert_eq!(entity_id_for_residue(&empty_subchain, &entities)?, ".");

        // Miss falls back to string_or_dot over the residue's label entity id.
        let miss = residue("ALA", EntityKind::Polymer, None, Some("Z"), Some("ownEnt"));
        assert_eq!(entity_id_for_residue(&miss, &entities)?, "ownEnt");
        let miss_empty = residue("ALA", EntityKind::Polymer, None, Some("Z"), Some(""));
        assert_eq!(entity_id_for_residue(&miss_empty, &entities)?, ".");
        let miss_absent = residue("ALA", EntityKind::Polymer, None, Some("Z"), None);
        assert_eq!(entity_id_for_residue(&miss_absent, &entities)?, ".");

        // Quoting is applied to both the entity name and the fallback.
        let spaced = entity("two words", EntityKind::Polymer, &["A"]);
        let in_a2 = residue("ALA", EntityKind::Polymer, None, Some("A"), None);
        assert_eq!(entity_id_for_residue(&in_a2, &[spaced])?, "'two words'");
        let spaced_fallback = residue("ALA", EntityKind::Polymer, None, None, Some("two words"));
        assert_eq!(entity_id_for_residue(&spaced_fallback, &[])?, "'two words'");
        Ok(())
    }

    // Silence unused warnings for construction helpers only used above.
    #[allow(dead_code)]
    fn _unused(parts: BioStructureParts) {
        let _ = parts;
        let _ = ChainKind::Protein;
        let _ = ChainSourceIds::default();
    }
}
