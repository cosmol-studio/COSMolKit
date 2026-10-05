//! Full selected-category mmCIF emission over canonical detached BIO data.
use super::value::{
    int_or_dot, int_or_qmark, number_or_dot, number_or_qmark, string_or_dot, string_or_qmark,
};
use super::{BioMmcifWriteError, BioMmcifWriteParams};
use crate::cif::CifBlock;
use cosmolkit_bio::*;

fn entity_kind_text(kind: EntityKind) -> &'static str {
    // Gemmi❗✔️: inline const char* entity_type_to_string(EntityType entity_type) {
    // Gemmi❗✔️:   switch (entity_type) {
    // Gemmi❗✔️:     case EntityType::Polymer: return "polymer";
    // Gemmi❗✔️:     case EntityType::Branched: return "branched";
    // Gemmi❗✔️:     case EntityType::NonPolymer: return "non-polymer";
    // Gemmi❗✔️:     case EntityType::Water: return "water";
    // Gemmi❗✔️:     default /*EntityType::Unknown*/: return "?";
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    match kind {
        EntityKind::Polymer => "polymer",
        EntityKind::Branched => "branched",
        EntityKind::NonPolymer => "non-polymer",
        EntityKind::Water => "water",
        EntityKind::Unknown => "?",
    }
}

fn polymer_kind_text(kind: PolymerKind) -> Option<&'static str> {
    // Gemmi❗✔️: inline const char* polymer_type_to_string(PolymerType polymer_type) {
    // Gemmi❗✔️:   switch (polymer_type) {
    // Gemmi❗✔️:     case PolymerType::PeptideL: return "polypeptide(L)";
    // Gemmi❗✔️:     case PolymerType::PeptideD: return "polypeptide(D)";
    // Gemmi❗✔️:     case PolymerType::Dna: return "polydeoxyribonucleotide";
    // Gemmi❗✔️:     case PolymerType::Rna: return "polyribonucleotide";
    // Gemmi❗✔️:     case PolymerType::DnaRnaHybrid:
    // Gemmi❗✔️:       return "'polydeoxyribonucleotide/polyribonucleotide hybrid'";
    // Gemmi❗✔️:     case PolymerType::SaccharideD: return "polysaccharide(D)";
    // Gemmi❗✔️:     case PolymerType::SaccharideL: return "polysaccharide(L)";
    // Gemmi❗✔️:     case PolymerType::Other: return "other";
    // Gemmi❗✔️:     case PolymerType::Pna: return "'peptide nucleic acid'";
    // Gemmi❗✔️:     case PolymerType::CyclicPseudoPeptide: return "cyclic-pseudo-peptide";
    // Gemmi❗✔️:     default /*PolymerType::Unknown*/: return "?";
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    match kind {
        PolymerKind::Unknown => None,
        PolymerKind::PeptideL => Some("polypeptide(L)"),
        PolymerKind::PeptideD => Some("polypeptide(D)"),
        PolymerKind::Dna => Some("polydeoxyribonucleotide"),
        PolymerKind::Rna => Some("polyribonucleotide"),
        PolymerKind::DnaRnaHybrid => Some("'polydeoxyribonucleotide/polyribonucleotide hybrid'"),
        PolymerKind::SaccharideD => Some("polysaccharide(D)"),
        PolymerKind::SaccharideL => Some("polysaccharide(L)"),
        PolymerKind::Pna => Some("'peptide nucleic acid'"),
        PolymerKind::CyclicPseudoPeptide => Some("cyclic-pseudo-peptide"),
        PolymerKind::Other => Some("other"),
    }
}

fn first_monomer(value: &str) -> &str {
    // Gemmi❗✔️:   static std::string first_mon(const std::string& mon_list) {
    // Gemmi❗✔️:     return mon_list.substr(0, mon_list.find(','));
    // Gemmi❗✔️:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    value.split_once(',').map_or(value, |(first, _)| first)
}

fn pdbx_one_letter_code(sequence: &[String], polymer_kind: PolymerKind) -> String {
    // Gemmi❗✔️: inline ResidueKind sequence_kind(PolymerType ptype) {
    // Gemmi❗✔️:   if (is_polypeptide(ptype))
    // Gemmi❗✔️:     return ResidueKind::AA;
    // Gemmi❗✔️:   if (ptype == PolymerType::Dna)
    // Gemmi❗✔️:     return ResidueKind::DNA;
    // Gemmi❗✔️:   if (ptype == PolymerType::Rna || ptype == PolymerType::DnaRnaHybrid)
    // Gemmi❗✔️:     return ResidueKind::RNA;
    // Gemmi❗✔️:   if (ptype == PolymerType::Unknown)
    // Gemmi❗✔️:     fail("sequence_kind(): unknown polymer type");
    // Gemmi❗✔️:   return ResidueKind::AA;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline bool is_polypeptide(PolymerType pt) {
    // Gemmi❗✔️:   return pt == PolymerType::PeptideL || pt == PolymerType::PeptideD;
    // Gemmi❗✔️: }
    // Behavior: the caller skips Unknown polymers before this private helper,
    // preserving the source dispatch guard before sequence_kind; known polymer
    // kinds use the copied DNA/RNA/polypeptide/default-AA mapping.
    // Complexity: this enum comparison is O(1), without allocation or row scans;
    // the existing string-building cost is reviewed separately below.
    // Gemmi❗❌: inline std::string pdbx_one_letter_code(const std::vector<std::string>& seq,
    // Gemmi❗❌:                                         ResidueKind kind) {
    // Gemmi❗❌:   std::string r;
    // Gemmi❗❌:   for (const std::string& item : seq) {
    // Gemmi❗❌:     std::string code = Entity::first_mon(item);
    // Gemmi❗❌:     const ResidueInfo ri = find_tabulated_residue(code);
    // Gemmi❗❌:     if (ri.is_standard() && ri.kind == kind)
    // Gemmi❗❌:       r += ri.one_letter_code;
    // Gemmi❗❌:     else
    // Gemmi❗❌:       cat_to(r, '(', code, ')');
    // Gemmi❗❌:   }
    // Gemmi❗❌:   return r;
    // Gemmi❗❌: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let mut result = String::new();
    for item in sequence {
        let code = first_monomer(item);
        let info = cosmolkit_bio::find_residue_info(code);
        let same_kind = info.kind
            == match polymer_kind {
                PolymerKind::Dna => ResidueInfoKind::Dna,
                PolymerKind::Rna | PolymerKind::DnaRnaHybrid => ResidueInfoKind::Rna,
                PolymerKind::PeptideL | PolymerKind::PeptideD => ResidueInfoKind::Aa,
                _ => ResidueInfoKind::Aa,
            };
        if info.is_standard() && same_kind {
            result.push(info.one_letter_code);
        } else {
            result.push('(');
            result.push_str(code);
            result.push(')');
        }
    }
    result
}

fn write_experiment_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
    entry_id: &str,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.exptl) {
    // Gemmi❗❌:     // _exptl
    // Gemmi❗❌:     if (!st.meta.experiments.empty()) {
    // Gemmi❗❌:       cif::Loop& loop = block.init_mmcif_loop("_exptl.",
    // Gemmi❗❌:                                               {"entry_id", "method", "crystals_number"});
    // Gemmi❗❌:       for (const ExperimentInfo& exper : st.meta.experiments)
    // Gemmi❗❌:         loop.add_row({id, cif::quote(exper.method),
    // Gemmi❗❌:                       int_or_qmark(exper.number_of_crystals)});
    // Gemmi❗❌:     } else {
    // Gemmi❗❌:       auto exptl_method = st.info.find("_exptl.method");
    // Gemmi❗❌:       if (exptl_method != st.info.end()) {
    // Gemmi❗❌:         cif::Loop& loop = block.init_mmcif_loop("_exptl.", {"entry_id", "method"});
    // Gemmi❗❌:         for (const std::string& m : gemmi::split_str(exptl_method->second, "; "))
    // Gemmi❗❌:           loop.add_row({id, cif::quote(m)});
    // Gemmi❗❌:       }
    // Gemmi❗❌:     }
    // Gemmi❗❌:
    // Gemmi❗❌:     // _exptl_crystal
    // Gemmi❗❌:     if (!st.meta.crystals.empty()) {
    // Gemmi❗❌:       cif::Loop& loop = block.init_mmcif_loop("_exptl_crystal.",
    // Gemmi❗❌:                                               {"id", "description"});
    // Gemmi❗❌:       for (const CrystalInfo& cryst : st.meta.crystals)
    // Gemmi❗❌:         loop.add_row({cryst.id, string_or_qmark(cryst.description)});
    // Gemmi❗❌:     }
    // Gemmi❗❌:
    // Gemmi❗❌:     // _exptl_crystal_grow
    // Gemmi❗❌:     if (std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
    // Gemmi❗❌:           [](const CrystalInfo& c) { return !c.ph_range.empty() || !std::isnan(c.ph); })) {
    // Gemmi❗❌:       cif::Loop& grow_loop = block.init_mmcif_loop("_exptl_crystal_grow.",
    // Gemmi❗❌:                                                    {"crystal_id", "pH", "pdbx_pH_range"});
    // Gemmi❗❌:       for (const CrystalInfo& crystal : st.meta.crystals)
    // Gemmi❗❌:         grow_loop.add_row({cif::quote(crystal.id),
    // Gemmi❗❌:                            number_or_qmark(crystal.ph),
    // Gemmi❗❌:                            string_or_qmark(crystal.ph_range)});
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if groups.exptl {
        if !structure.metadata().experiments.is_empty() {
            let rows = structure
                .metadata()
                .experiments
                .iter()
                .map(|experiment| {
                    vec![
                        entry_id.to_string(),
                        quote_cif_value(&experiment.method),
                        int_or_qmark(Some(experiment.number_of_crystals)),
                    ]
                })
                .collect();
            add_mmcif_rows(
                block,
                "_exptl.",
                &["entry_id", "method", "crystals_number"],
                rows,
            )?;
        } else if let Some(methods) = structure.source_state().info.get("_exptl.method") {
            let rows = methods
                .split("; ")
                .map(|method| vec![entry_id.to_string(), quote_cif_value(method)])
                .collect();
            add_mmcif_rows(block, "_exptl.", &["entry_id", "method"], rows)?;
        }
        if !structure.metadata().crystals.is_empty() {
            let rows = structure
                .metadata()
                .crystals
                .iter()
                .map(|crystal| vec![crystal.id.clone(), string_or_qmark(&crystal.description)])
                .collect();
            add_mmcif_rows(block, "_exptl_crystal.", &["id", "description"], rows)?;
        }
        if structure
            .metadata()
            .crystals
            .iter()
            .any(|crystal| !crystal.ph_range.is_empty() || !crystal.ph.is_nan())
        {
            let rows = structure
                .metadata()
                .crystals
                .iter()
                .map(|crystal| {
                    vec![
                        quote_cif_value(&crystal.id),
                        number_or_qmark(Some(crystal.ph)),
                        string_or_qmark(&crystal.ph_range),
                    ]
                })
                .collect();
            add_mmcif_rows(
                block,
                "_exptl_crystal_grow.",
                &["crystal_id", "pH", "pdbx_pH_range"],
                rows,
            )?;
        }
    }
    if groups.diffrn
        && structure
            .metadata()
            .crystals
            .iter()
            .any(|crystal| !crystal.diffractions.is_empty())
    {
        write_diffraction_categories(structure, block)?;
    }
    Ok(())
}

fn write_diffraction_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.diffrn &&
    // Gemmi❗❌:       std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
    // Gemmi❗❌:                   [](const CrystalInfo& c) { return !c.diffractions.empty(); })) {
    // Gemmi❗❌:
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_diffrn.", {"id", "crystal_id", "ambient_temp"});
    // Gemmi❗❌:     for (const CrystalInfo& cryst : st.meta.crystals)
    // Gemmi❗❌:       for (const DiffractionInfo& diffr : cryst.diffractions)
    // Gemmi❗❌:         loop.add_row({diffr.id, cryst.id, number_or_qmark(diffr.temperature)});
    // Gemmi❗❌:     // _diffrn_detector
    // Gemmi❗❌:     cif::Loop& det_loop = block.init_mmcif_loop("_diffrn_detector.",
    // Gemmi❗❌:                                                 {"diffrn_id",
    // Gemmi❗❌:                                                  "pdbx_collection_date",
    // Gemmi❗❌:                                                  "detector",
    // Gemmi❗❌:                                                  "type",
    // Gemmi❗❌:                                                  "details"});
    // Gemmi❗❌:     for (const CrystalInfo& cryst : st.meta.crystals)
    // Gemmi❗❌:       for (const DiffractionInfo& diffr : cryst.diffractions)
    // Gemmi❗❌:         det_loop.add_row({diffr.id,
    // Gemmi❗❌:                           string_or_qmark(diffr.collection_date),
    // Gemmi❗❌:                           string_or_qmark(diffr.detector),
    // Gemmi❗❌:                           string_or_qmark(diffr.detector_make),
    // Gemmi❗❌:                           string_or_qmark(diffr.optics)});
    // Gemmi❗❌:
    // Gemmi❗❌:     // _diffrn_radiation
    // Gemmi❗❌:     cif::Loop& rad_loop = block.init_mmcif_loop("_diffrn_radiation.",
    // Gemmi❗❌:                                                 {"diffrn_id",
    // Gemmi❗❌:                                                  "pdbx_scattering_type",
    // Gemmi❗❌:                                                  "pdbx_monochromatic_or_laue_m_l",
    // Gemmi❗❌:                                                  "monochromator"});
    // Gemmi❗❌:     for (const CrystalInfo& cryst : st.meta.crystals)
    // Gemmi❗❌:       for (const DiffractionInfo& diffr : cryst.diffractions)
    // Gemmi❗❌:         rad_loop.add_row({diffr.id,
    // Gemmi❗❌:                           string_or_qmark(diffr.scattering_type),
    // Gemmi❗❌:                           std::string(1, diffr.mono_or_laue ? diffr.mono_or_laue : '?'),
    // Gemmi❗❌:                           string_or_qmark(diffr.monochromator)});
    // Gemmi❗❌:     // _diffrn_source
    // Gemmi❗❌:     cif::Loop& source_loop = block.init_mmcif_loop("_diffrn_source.",
    // Gemmi❗❌:                                                    {"diffrn_id",
    // Gemmi❗❌:                                                     "source",
    // Gemmi❗❌:                                                     "type",
    // Gemmi❗❌:                                                     "pdbx_synchrotron_site",
    // Gemmi❗❌:                                                     "pdbx_synchrotron_beamline",
    // Gemmi❗❌:                                                     "pdbx_wavelength_list"});
    // Gemmi❗❌:     for (const CrystalInfo& crystal : st.meta.crystals)
    // Gemmi❗❌:       for (const DiffractionInfo& diffr : crystal.diffractions)
    // Gemmi❗❌:         source_loop.add_row({diffr.id,
    // Gemmi❗❌:                              string_or_qmark(diffr.source),
    // Gemmi❗❌:                              string_or_qmark(diffr.source_type),
    // Gemmi❗❌:                              string_or_qmark(diffr.synchrotron),
    // Gemmi❗❌:                              string_or_qmark(diffr.beamline),
    // Gemmi❗❌:                              string_or_qmark(diffr.wavelengths)});
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let diffractions = structure
        .metadata()
        .crystals
        .iter()
        .flat_map(|crystal| {
            crystal
                .diffractions
                .iter()
                .map(move |diffraction| (crystal, diffraction))
        })
        .collect::<Vec<_>>();
    add_mmcif_rows(
        block,
        "_diffrn.",
        &["id", "crystal_id", "ambient_temp"],
        diffractions
            .iter()
            .map(|(crystal, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    crystal.id.clone(),
                    number_or_qmark(Some(diffraction.temperature)),
                ]
            })
            .collect(),
    )?;
    add_mmcif_rows(
        block,
        "_diffrn_detector.",
        &[
            "diffrn_id",
            "pdbx_collection_date",
            "detector",
            "type",
            "details",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.collection_date),
                    string_or_qmark(&diffraction.detector),
                    string_or_qmark(&diffraction.detector_make),
                    string_or_qmark(&diffraction.optics),
                ]
            })
            .collect(),
    )?;
    add_mmcif_rows(
        block,
        "_diffrn_radiation.",
        &[
            "diffrn_id",
            "pdbx_scattering_type",
            "pdbx_monochromatic_or_laue_m_l",
            "monochromator",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.scattering_type),
                    char::from(if diffraction.mono_or_laue == 0 {
                        b'?'
                    } else {
                        diffraction.mono_or_laue
                    })
                    .to_string(),
                    string_or_qmark(&diffraction.monochromator),
                ]
            })
            .collect(),
    )?;
    add_mmcif_rows(
        block,
        "_diffrn_source.",
        &[
            "diffrn_id",
            "source",
            "type",
            "pdbx_synchrotron_site",
            "pdbx_synchrotron_beamline",
            "pdbx_wavelength_list",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.source),
                    string_or_qmark(&diffraction.source_type),
                    string_or_qmark(&diffraction.synchrotron),
                    string_or_qmark(&diffraction.beamline),
                    string_or_qmark(&diffraction.wavelengths),
                ]
            })
            .collect(),
    )
}

fn write_reflection_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
    entry_id: &str,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.reflns && !st.meta.experiments.empty()) {
    // Gemmi❗❌:     // _reflns
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_reflns.", {
    // Gemmi❗❌:         "entry_id",
    // Gemmi❗❌:         "pdbx_ordinal",
    // Gemmi❗❌:         "pdbx_diffrn_id",
    // Gemmi❗❌:         "number_obs",
    // Gemmi❗❌:         "d_resolution_high",
    // Gemmi❗❌:         "d_resolution_low",
    // Gemmi❗❌:         "percent_possible_obs",
    // Gemmi❗❌:         "pdbx_redundancy",
    // Gemmi❗❌:         "pdbx_Rmerge_I_obs",
    // Gemmi❗❌:         "pdbx_Rsym_value",
    // Gemmi❗❌:         "pdbx_netI_over_sigmaI",
    // Gemmi❗❌:         /*"B_iso_Wilson_estimate"*/});
    // Gemmi❗❌:     int n = 0;
    // Gemmi❗❌:     for (const ExperimentInfo& exper : st.meta.experiments)
    // Gemmi❗❌:       loop.add_row({id,
    // Gemmi❗❌:                     std::to_string(++n),
    // Gemmi❗❌:                     string_or_dot(join_str(exper.diffraction_ids, ",")),
    // Gemmi❗❌:                     int_or_qmark(exper.unique_reflections),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.resolution_high),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.resolution_low),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.completeness),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.redundancy),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.r_merge),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.r_sym),
    // Gemmi❗❌:                     number_or_qmark(exper.reflections.mean_I_over_sigma),
    // Gemmi❗❌:                     /*number_or_qmark(exper.b_wilson)*/});
    // Gemmi❗❌:     // _reflns_shell
    // Gemmi❗❌:     cif::Loop* shell_loop = nullptr;
    // Gemmi❗❌:     n = 0;
    // Gemmi❗❌:     for (const ExperimentInfo& exper : st.meta.experiments) {
    // Gemmi❗❌:       std::string diffrn_id =
    // Gemmi❗❌:         string_or_dot(join_str(exper.diffraction_ids, ","));
    // Gemmi❗❌:       for (const ReflectionsInfo& shell : exper.shells) {
    // Gemmi❗❌:         if (!shell_loop)
    // Gemmi❗❌:           shell_loop = &block.init_mmcif_loop("_reflns_shell.", {
    // Gemmi❗❌:               "pdbx_ordinal",
    // Gemmi❗❌:               "pdbx_diffrn_id",
    // Gemmi❗❌:               "d_res_high",
    // Gemmi❗❌:               "d_res_low",
    // Gemmi❗❌:               "percent_possible_all",
    // Gemmi❗❌:               "pdbx_redundancy",
    // Gemmi❗❌:               "Rmerge_I_obs",
    // Gemmi❗❌:               "pdbx_Rsym_value",
    // Gemmi❗❌:               "meanI_over_sigI_obs"});
    // Gemmi❗❌:
    // Gemmi❗❌:         shell_loop->add_row({std::to_string(++n),
    // Gemmi❗❌:                              diffrn_id,
    // Gemmi❗❌:                              number_or_qmark(shell.resolution_high),
    // Gemmi❗❌:                              number_or_qmark(shell.resolution_low),
    // Gemmi❗❌:                              number_or_qmark(shell.completeness),
    // Gemmi❗❌:                              number_or_qmark(shell.redundancy),
    // Gemmi❗❌:                              number_or_qmark(shell.r_merge),
    // Gemmi❗❌:                              number_or_qmark(shell.r_sym),
    // Gemmi❗❌:                              number_or_qmark(shell.mean_I_over_sigma)});
    // Gemmi❗❌:       }
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if !groups.reflns || structure.metadata().experiments.is_empty() {
        return Ok(());
    }
    let rows = structure
        .metadata()
        .experiments
        .iter()
        .enumerate()
        .map(|(index, experiment)| {
            vec![
                entry_id.to_string(),
                (index + 1).to_string(),
                string_or_dot(&experiment.diffraction_ids.join(",")),
                int_or_qmark(Some(experiment.unique_reflections)),
                number_or_qmark(Some(experiment.reflections.resolution_high)),
                number_or_qmark(Some(experiment.reflections.resolution_low)),
                number_or_qmark(Some(experiment.reflections.completeness)),
                number_or_qmark(Some(experiment.reflections.redundancy)),
                number_or_qmark(Some(experiment.reflections.r_merge)),
                number_or_qmark(Some(experiment.reflections.r_sym)),
                number_or_qmark(Some(experiment.reflections.mean_i_over_sigma)),
            ]
        })
        .collect();
    add_mmcif_rows(
        block,
        "_reflns.",
        &[
            "entry_id",
            "pdbx_ordinal",
            "pdbx_diffrn_id",
            "number_obs",
            "d_resolution_high",
            "d_resolution_low",
            "percent_possible_obs",
            "pdbx_redundancy",
            "pdbx_Rmerge_I_obs",
            "pdbx_Rsym_value",
            "pdbx_netI_over_sigmaI",
        ],
        rows,
    )?;
    let mut shell_rows = Vec::new();
    let mut ordinal = 0_usize;
    for experiment in &structure.metadata().experiments {
        let diffraction_id = string_or_dot(&experiment.diffraction_ids.join(","));
        for shell in &experiment.shells {
            ordinal += 1;
            shell_rows.push(vec![
                ordinal.to_string(),
                diffraction_id.clone(),
                number_or_qmark(Some(shell.resolution_high)),
                number_or_qmark(Some(shell.resolution_low)),
                number_or_qmark(Some(shell.completeness)),
                number_or_qmark(Some(shell.redundancy)),
                number_or_qmark(Some(shell.r_merge)),
                number_or_qmark(Some(shell.r_sym)),
                number_or_qmark(Some(shell.mean_i_over_sigma)),
            ]);
        }
    }
    if !shell_rows.is_empty() {
        add_mmcif_rows(
            block,
            "_reflns_shell.",
            &[
                "pdbx_ordinal",
                "pdbx_diffrn_id",
                "d_res_high",
                "d_res_low",
                "percent_possible_all",
                "pdbx_redundancy",
                "Rmerge_I_obs",
                "pdbx_Rsym_value",
                "meanI_over_sigI_obs",
            ],
            shell_rows,
        )?;
    }
    Ok(())
}

fn write_refinement_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
    entry_id: &str,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.refine && !st.meta.refinement.empty()) {
    // Gemmi❗❌:     block.items.reserve(block.items.size() + 4);
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_refine.", {
    // Gemmi❗❌:         "entry_id",
    // Gemmi❗❌:         "pdbx_refine_id",
    // Gemmi❗❌:         "ls_d_res_high",
    // Gemmi❗❌:         "ls_d_res_low",
    // Gemmi❗❌:         "ls_percent_reflns_obs",
    // Gemmi❗❌:         "ls_number_reflns_obs",
    // Gemmi❗❌:         "ls_number_reflns_R_work"});
    // Gemmi❗❌:     cif::Loop& analyze_loop = block.init_mmcif_loop("_refine_analyze.", {
    // Gemmi❗❌:         "entry_id",
    // Gemmi❗❌:         "pdbx_refine_id",
    // Gemmi❗❌:         "Luzzati_coordinate_error_obs"});
    // Gemmi❗❌:     cif::Loop& restr_loop = block.init_mmcif_loop("_refine_ls_restr.", {
    // Gemmi❗❌:         "pdbx_refine_id", "type",
    // Gemmi❗❌:         "number", "weight", "pdbx_restraint_function", "dev_ideal"});
    // Gemmi❗❌:     // _refine_ls_shell
    // Gemmi❗❌:     std::vector<std::string> shell_tags = {
    // Gemmi❗❌:         "pdbx_refine_id",
    // Gemmi❗❌:         "d_res_high",
    // Gemmi❗❌:         "d_res_low",
    // Gemmi❗❌:         "percent_reflns_obs",
    // Gemmi❗❌:         "number_reflns_obs",
    // Gemmi❗❌:         "number_reflns_R_work",
    // Gemmi❗❌:         "number_reflns_R_free",
    // Gemmi❗❌:         "R_factor_obs",
    // Gemmi❗❌:         "R_factor_R_work",
    // Gemmi❗❌:         "R_factor_R_free"};
    // Gemmi❗❌:     bool has_shell_fsc = false;
    // Gemmi❗❌:     bool has_shell_ffcc = false;
    // Gemmi❗❌:     bool has_shell_iicc = false;
    // Gemmi❗❌:     for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi❗❌:       for (const BasicRefinementInfo& bin : ref.bins) {
    // Gemmi❗❌:         if (!std::isnan(bin.fsc_work) || !std::isnan(bin.fsc_free))
    // Gemmi❗❌:           has_shell_fsc = true;
    // Gemmi❗❌:         if (!std::isnan(bin.cc_fo_fc_work) || !std::isnan(bin.cc_fo_fc_free))
    // Gemmi❗❌:           has_shell_ffcc = true;
    // Gemmi❗❌:         if (!std::isnan(bin.cc_intensity_work) || !std::isnan(bin.cc_intensity_free))
    // Gemmi❗❌:           has_shell_iicc = true;
    // Gemmi❗❌:       }
    // Gemmi❗❌:     if (has_shell_fsc) {
    // Gemmi❗❌:       shell_tags.push_back("pdbx_fsc_work");
    // Gemmi❗❌:       shell_tags.push_back("pdbx_fsc_free");
    // Gemmi❗❌:     }
    // Gemmi❗❌:     if (has_shell_ffcc) {
    // Gemmi❗❌:       shell_tags.push_back("correlation_coeff_Fo_to_Fc");
    // Gemmi❗❌:       shell_tags.push_back("correlation_coeff_Fo_to_Fc_free");
    // Gemmi❗❌:     }
    // Gemmi❗❌:     if (has_shell_iicc) {
    // Gemmi❗❌:       shell_tags.push_back("correlation_coeff_I_to_Fcsqd_work");
    // Gemmi❗❌:       shell_tags.push_back("correlation_coeff_I_to_Fcsqd_free");
    // Gemmi❗❌:     }
    // Gemmi❗❌:     cif::Loop& shell_loop = block.init_mmcif_loop("_refine_ls_shell.", shell_tags);
    // Gemmi❗❌:
    // Gemmi❗❌:     for (size_t i = 0; i != st.meta.refinement.size(); ++i) {
    // Gemmi❗❌:       const RefinementInfo& ref = st.meta.refinement[i];
    // Gemmi❗❌:       loop.add_values({id,
    // Gemmi❗❌:                        cif::quote(ref.id),
    // Gemmi❗❌:                        number_or_dot(ref.resolution_high),
    // Gemmi❗❌:                        number_or_dot(ref.resolution_low),
    // Gemmi❗❌:                        number_or_dot(ref.completeness),
    // Gemmi❗❌:                        int_or_dot(get_number_obs(ref)),
    // Gemmi❗❌:                        int_or_qmark(get_number_work(ref))});
    // Gemmi❗❌:       auto add = [&](const std::string& tag, const std::string& val) {
    // Gemmi❗❌:         if (i == 0)
    // Gemmi❗❌:           loop.tags.push_back("_refine." + tag);
    // Gemmi❗❌:         loop.values.push_back(val);
    // Gemmi❗❌:       };
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::rfree_set_count))
    // Gemmi❗❌:         add("ls_number_reflns_R_free", int_or_dot(ref.rfree_set_count));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::r_all))
    // Gemmi❗❌:         add("ls_R_factor_obs", number_or_qmark(ref.r_all));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::r_work))
    // Gemmi❗❌:         add("ls_R_factor_R_work", number_or_qmark(ref.r_work));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::r_free))
    // Gemmi❗❌:         add("ls_R_factor_R_free", number_or_qmark(ref.r_free));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::cross_validation_method))
    // Gemmi❗❌:         add("pdbx_ls_cross_valid_method",
    // Gemmi❗❌:             string_or_qmark(ref.cross_validation_method));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::rfree_selection_method))
    // Gemmi❗❌:         add("pdbx_R_Free_selection_details",
    // Gemmi❗❌:             string_or_qmark(ref.rfree_selection_method));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::mean_b))
    // Gemmi❗❌:         add("B_iso_mean", number_or_qmark(ref.mean_b));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::aniso_b)) {
    // Gemmi❗❌:         if (i == 0)
    // Gemmi❗❌:           for (const char* index : {"[1][1]", "[2][2]", "[3][3]", "[1][2]", "[1][3]", "[2][3]"})
    // Gemmi❗❌:             loop.tags.push_back(std::string("_refine.aniso_B") + index);
    // Gemmi❗❌:         for (double d : ref.aniso_b.elements_pdb())
    // Gemmi❗❌:           loop.values.push_back(number_or_qmark(d));
    // Gemmi❗❌:       }
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::dpi_blow_r))
    // Gemmi❗❌:         add("pdbx_overall_SU_R_Blow_DPI",
    // Gemmi❗❌:             number_or_qmark(ref.dpi_blow_r));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::dpi_blow_rfree))
    // Gemmi❗❌:         add("pdbx_overall_SU_R_free_Blow_DPI",
    // Gemmi❗❌:             number_or_qmark(ref.dpi_blow_rfree));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::dpi_cruickshank_r))
    // Gemmi❗❌:         add("overall_SU_R_Cruickshank_DPI",
    // Gemmi❗❌:             number_or_qmark(ref.dpi_cruickshank_r));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::dpi_cruickshank_rfree))
    // Gemmi❗❌:         add("pdbx_overall_SU_R_free_Cruickshank_DPI",
    // Gemmi❗❌:             number_or_qmark(ref.dpi_cruickshank_rfree));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::cc_fo_fc_work))
    // Gemmi❗❌:         add("correlation_coeff_Fo_to_Fc", number_or_qmark(ref.cc_fo_fc_work));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::cc_fo_fc_free))
    // Gemmi❗❌:         add("correlation_coeff_Fo_to_Fc_free",
    // Gemmi❗❌:             number_or_qmark(ref.cc_fo_fc_free));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::fsc_work))
    // Gemmi❗❌:         add("pdbx_average_fsc_work", number_or_qmark(ref.fsc_work));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::fsc_free))
    // Gemmi❗❌:         add("pdbx_average_fsc_free", number_or_qmark(ref.fsc_free));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::cc_intensity_work))
    // Gemmi❗❌:         add("correlation_coeff_I_to_Fcsqd_work", number_or_qmark(ref.cc_intensity_work));
    // Gemmi❗❌:       if (st.meta.has(&RefinementInfo::cc_intensity_free))
    // Gemmi❗❌:         add("correlation_coeff_I_to_Fcsqd_free", number_or_qmark(ref.cc_intensity_free));
    // Gemmi❗❌:       if (!st.meta.solved_by.empty())
    // Gemmi❗❌:         add("pdbx_method_to_determine_struct", string_or_qmark(st.meta.solved_by));
    // Gemmi❗❌:       if (!st.meta.starting_model.empty())
    // Gemmi❗❌:         add("pdbx_starting_model", string_or_qmark(st.meta.starting_model));
    // Gemmi❗❌:       if (!std::isnan(ref.luzzati_error))
    // Gemmi❗❌:         analyze_loop.add_row({id,
    // Gemmi❗❌:                               cif::quote(ref.id),
    // Gemmi❗❌:                               number_or_qmark(ref.luzzati_error)});
    // Gemmi❗❌:       for (const RefinementInfo::Restr& restr : ref.restr_stats)
    // Gemmi❗❌:         restr_loop.add_row({cif::quote(ref.id),
    // Gemmi❗❌:                             cif::quote(restr.name),
    // Gemmi❗❌:                             int_or_qmark(restr.count),
    // Gemmi❗❌:                             number_or_qmark(restr.weight),
    // Gemmi❗❌:                             string_or_qmark(restr.function),
    // Gemmi❗❌:                             number_or_qmark(restr.dev_ideal)});
    // Gemmi❗❌:       for (const BasicRefinementInfo& bin : ref.bins) {
    // Gemmi❗❌:         shell_loop.add_values({cif::quote(ref.id),
    // Gemmi❗❌:                                number_or_dot(bin.resolution_high),
    // Gemmi❗❌:                                number_or_qmark(bin.resolution_low),
    // Gemmi❗❌:                                number_or_qmark(bin.completeness),
    // Gemmi❗❌:                                int_or_qmark(get_number_obs(bin)),
    // Gemmi❗❌:                                int_or_qmark(get_number_work(bin)),
    // Gemmi❗❌:                                int_or_qmark(bin.rfree_set_count),
    // Gemmi❗❌:                                number_or_qmark(bin.r_all),
    // Gemmi❗❌:                                number_or_qmark(bin.r_work),
    // Gemmi❗❌:                                number_or_qmark(bin.r_free)});
    // Gemmi❗❌:         if (has_shell_fsc)
    // Gemmi❗❌:           shell_loop.add_values({number_or_qmark(bin.fsc_work),
    // Gemmi❗❌:                                  number_or_qmark(bin.fsc_free)});
    // Gemmi❗❌:         if (has_shell_ffcc)
    // Gemmi❗❌:           shell_loop.add_values({number_or_qmark(bin.cc_fo_fc_work),
    // Gemmi❗❌:                                  number_or_qmark(bin.cc_fo_fc_free)});
    // Gemmi❗❌:         if (has_shell_iicc)
    // Gemmi❗❌:           shell_loop.add_values({number_or_qmark(bin.cc_intensity_work),
    // Gemmi❗❌:                                  number_or_qmark(bin.cc_intensity_free)});
    // Gemmi❗❌:       }
    // Gemmi❗❌:     }
    // Gemmi❗❌:     assert(shell_loop.values.size() % shell_loop.tags.size() == 0);
    // Gemmi❗❌:     assert(loop.values.size() % loop.tags.size() == 0);
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if !groups.refine || structure.metadata().refinement.is_empty() {
        return Ok(());
    }
    let bins = structure
        .metadata()
        .refinement
        .iter()
        .flat_map(|refinement| &refinement.bins)
        .collect::<Vec<_>>();
    let has_shell_fsc = bins
        .iter()
        .any(|bin| !bin.fsc_work.is_nan() || !bin.fsc_free.is_nan());
    let has_shell_ffcc = bins
        .iter()
        .any(|bin| !bin.cc_fo_fc_work.is_nan() || !bin.cc_fo_fc_free.is_nan());
    let has_shell_iicc = bins
        .iter()
        .any(|bin| !bin.cc_intensity_work.is_nan() || !bin.cc_intensity_free.is_nan());
    let refinements = &structure.metadata().refinement;
    let has_rfree_count = refinements
        .iter()
        .any(|value| (value.basic.rfree_set_count != -1));
    let has_r_all = refinements.iter().any(|value| !value.basic.r_all.is_nan());
    let has_r_work = refinements.iter().any(|value| !value.basic.r_work.is_nan());
    let has_r_free = refinements.iter().any(|value| !value.basic.r_free.is_nan());
    let has_cross_validation = refinements
        .iter()
        .any(|value| !value.cross_validation_method.is_empty());
    let has_rfree_selection = refinements
        .iter()
        .any(|value| !value.rfree_selection_method.is_empty());
    let has_mean_b = refinements.iter().any(|value| !value.mean_b.is_nan());
    let has_aniso = refinements.iter().any(|value| !value.aniso_b[0].is_nan());
    let has_dpi_blow_r = refinements.iter().any(|value| !value.dpi_blow_r.is_nan());
    let has_dpi_blow_rfree = refinements
        .iter()
        .any(|value| !value.dpi_blow_rfree.is_nan());
    let has_dpi_cruickshank_r = refinements
        .iter()
        .any(|value| !value.dpi_cruickshank_r.is_nan());
    let has_dpi_cruickshank_rfree = refinements
        .iter()
        .any(|value| !value.dpi_cruickshank_rfree.is_nan());
    let has_cc_fo_fc_work = refinements
        .iter()
        .any(|value| !value.basic.cc_fo_fc_work.is_nan());
    let has_cc_fo_fc_free = refinements
        .iter()
        .any(|value| !value.basic.cc_fo_fc_free.is_nan());
    let has_fsc_work = refinements
        .iter()
        .any(|value| !value.basic.fsc_work.is_nan());
    let has_fsc_free = refinements
        .iter()
        .any(|value| !value.basic.fsc_free.is_nan());
    let has_cc_intensity_work = refinements
        .iter()
        .any(|value| !value.basic.cc_intensity_work.is_nan());
    let has_cc_intensity_free = refinements
        .iter()
        .any(|value| !value.basic.cc_intensity_free.is_nan());

    let mut refine_tags = vec![
        "entry_id",
        "pdbx_refine_id",
        "ls_d_res_high",
        "ls_d_res_low",
        "ls_percent_reflns_obs",
        "ls_number_reflns_obs",
        "ls_number_reflns_R_work",
    ];
    let optional_tags = [
        (has_rfree_count, "ls_number_reflns_R_free"),
        (has_r_all, "ls_R_factor_obs"),
        (has_r_work, "ls_R_factor_R_work"),
        (has_r_free, "ls_R_factor_R_free"),
        (has_cross_validation, "pdbx_ls_cross_valid_method"),
        (has_rfree_selection, "pdbx_R_Free_selection_details"),
        (has_mean_b, "B_iso_mean"),
    ];
    for (present, tag) in optional_tags {
        if present {
            refine_tags.push(tag);
        }
    }
    if has_aniso {
        refine_tags.extend([
            "aniso_B[1][1]",
            "aniso_B[2][2]",
            "aniso_B[3][3]",
            "aniso_B[1][2]",
            "aniso_B[1][3]",
            "aniso_B[2][3]",
        ]);
    }
    for (present, tag) in [
        (has_dpi_blow_r, "pdbx_overall_SU_R_Blow_DPI"),
        (has_dpi_blow_rfree, "pdbx_overall_SU_R_free_Blow_DPI"),
        (has_dpi_cruickshank_r, "overall_SU_R_Cruickshank_DPI"),
        (
            has_dpi_cruickshank_rfree,
            "pdbx_overall_SU_R_free_Cruickshank_DPI",
        ),
        (has_cc_fo_fc_work, "correlation_coeff_Fo_to_Fc"),
        (has_cc_fo_fc_free, "correlation_coeff_Fo_to_Fc_free"),
        (has_fsc_work, "pdbx_average_fsc_work"),
        (has_fsc_free, "pdbx_average_fsc_free"),
        (has_cc_intensity_work, "correlation_coeff_I_to_Fcsqd_work"),
        (has_cc_intensity_free, "correlation_coeff_I_to_Fcsqd_free"),
    ] {
        if present {
            refine_tags.push(tag);
        }
    }
    let has_solved_by = !structure.metadata().solved_by.is_empty();
    let has_starting_model = !structure.metadata().starting_model.is_empty();
    if has_solved_by {
        refine_tags.push("pdbx_method_to_determine_struct");
    }
    if has_starting_model {
        refine_tags.push("pdbx_starting_model");
    }

    let mut refine_rows = Vec::with_capacity(refinements.len());
    let mut analyze_rows = Vec::new();
    let mut restraint_rows = Vec::new();
    let mut shell_rows = Vec::new();
    for refinement in refinements {
        let mut row = vec![
            entry_id.to_string(),
            quote_cif_value(&refinement.id),
            number_or_dot(Some(refinement.basic.resolution_high)),
            number_or_dot(Some(refinement.basic.resolution_low)),
            number_or_dot(Some(refinement.basic.completeness)),
            int_or_dot(get_number_obs(&refinement.basic)),
            int_or_qmark(get_number_work(&refinement.basic)),
        ];
        if has_rfree_count {
            row.push(int_or_dot(Some(refinement.basic.rfree_set_count)));
        }
        if has_r_all {
            row.push(number_or_qmark(Some(refinement.basic.r_all)));
        }
        if has_r_work {
            row.push(number_or_qmark(Some(refinement.basic.r_work)));
        }
        if has_r_free {
            row.push(number_or_qmark(Some(refinement.basic.r_free)));
        }
        if has_cross_validation {
            row.push(string_or_qmark(&refinement.cross_validation_method));
        }
        if has_rfree_selection {
            row.push(string_or_qmark(&refinement.rfree_selection_method));
        }
        if has_mean_b {
            row.push(number_or_qmark(Some(refinement.mean_b)));
        }
        if has_aniso {
            row.extend([
                number_or_qmark(Some(refinement.aniso_b[0])),
                number_or_qmark(Some(refinement.aniso_b[1])),
                number_or_qmark(Some(refinement.aniso_b[2])),
                number_or_qmark(Some(refinement.aniso_b[3])),
                number_or_qmark(Some(refinement.aniso_b[4])),
                number_or_qmark(Some(refinement.aniso_b[5])),
            ]);
        }
        if has_dpi_blow_r {
            row.push(number_or_qmark(Some(refinement.dpi_blow_r)));
        }
        if has_dpi_blow_rfree {
            row.push(number_or_qmark(Some(refinement.dpi_blow_rfree)));
        }
        if has_dpi_cruickshank_r {
            row.push(number_or_qmark(Some(refinement.dpi_cruickshank_r)));
        }
        if has_dpi_cruickshank_rfree {
            row.push(number_or_qmark(Some(refinement.dpi_cruickshank_rfree)));
        }
        if has_cc_fo_fc_work {
            row.push(number_or_qmark(Some(refinement.basic.cc_fo_fc_work)));
        }
        if has_cc_fo_fc_free {
            row.push(number_or_qmark(Some(refinement.basic.cc_fo_fc_free)));
        }
        if has_fsc_work {
            row.push(number_or_qmark(Some(refinement.basic.fsc_work)));
        }
        if has_fsc_free {
            row.push(number_or_qmark(Some(refinement.basic.fsc_free)));
        }
        if has_cc_intensity_work {
            row.push(number_or_qmark(Some(refinement.basic.cc_intensity_work)));
        }
        if has_cc_intensity_free {
            row.push(number_or_qmark(Some(refinement.basic.cc_intensity_free)));
        }
        if has_solved_by {
            row.push(string_or_qmark(&structure.metadata().solved_by));
        }
        if has_starting_model {
            row.push(string_or_qmark(&structure.metadata().starting_model));
        }
        refine_rows.push(row);

        if !refinement.luzzati_error.is_nan() {
            analyze_rows.push(vec![
                entry_id.to_string(),
                quote_cif_value(&refinement.id),
                number_or_qmark(Some(refinement.luzzati_error)),
            ]);
        }
        for restraint in &refinement.restr_stats {
            restraint_rows.push(vec![
                quote_cif_value(&refinement.id),
                quote_cif_value(&restraint.name),
                int_or_qmark(Some(restraint.count)),
                number_or_qmark(Some(restraint.weight)),
                string_or_qmark(&restraint.function),
                number_or_qmark(Some(restraint.dev_ideal)),
            ]);
        }
        for bin in &refinement.bins {
            let mut shell = vec![
                quote_cif_value(&refinement.id),
                number_or_dot(Some(bin.resolution_high)),
                number_or_qmark(Some(bin.resolution_low)),
                number_or_qmark(Some(bin.completeness)),
                int_or_qmark(get_number_obs(bin)),
                int_or_qmark(get_number_work(bin)),
                int_or_qmark(Some(bin.rfree_set_count)),
                number_or_qmark(Some(bin.r_all)),
                number_or_qmark(Some(bin.r_work)),
                number_or_qmark(Some(bin.r_free)),
            ];
            if has_shell_fsc {
                shell.extend([
                    number_or_qmark(Some(bin.fsc_work)),
                    number_or_qmark(Some(bin.fsc_free)),
                ]);
            }
            if has_shell_ffcc {
                shell.extend([
                    number_or_qmark(Some(bin.cc_fo_fc_work)),
                    number_or_qmark(Some(bin.cc_fo_fc_free)),
                ]);
            }
            if has_shell_iicc {
                shell.extend([
                    number_or_qmark(Some(bin.cc_intensity_work)),
                    number_or_qmark(Some(bin.cc_intensity_free)),
                ]);
            }
            shell_rows.push(shell);
        }
    }
    add_mmcif_rows(block, "_refine.", &refine_tags, refine_rows)?;
    add_mmcif_rows(
        block,
        "_refine_analyze.",
        &["entry_id", "pdbx_refine_id", "Luzzati_coordinate_error_obs"],
        analyze_rows,
    )?;
    add_mmcif_rows(
        block,
        "_refine_ls_restr.",
        &[
            "pdbx_refine_id",
            "type",
            "number",
            "weight",
            "pdbx_restraint_function",
            "dev_ideal",
        ],
        restraint_rows,
    )?;
    let mut shell_tags = vec![
        "pdbx_refine_id",
        "d_res_high",
        "d_res_low",
        "percent_reflns_obs",
        "number_reflns_obs",
        "number_reflns_R_work",
        "number_reflns_R_free",
        "R_factor_obs",
        "R_factor_R_work",
        "R_factor_R_free",
    ];
    if has_shell_fsc {
        shell_tags.extend(["pdbx_fsc_work", "pdbx_fsc_free"]);
    }
    if has_shell_ffcc {
        shell_tags.extend([
            "correlation_coeff_Fo_to_Fc",
            "correlation_coeff_Fo_to_Fc_free",
        ]);
    }
    if has_shell_iicc {
        shell_tags.extend([
            "correlation_coeff_I_to_Fcsqd_work",
            "correlation_coeff_I_to_Fcsqd_free",
        ]);
    }
    add_mmcif_rows(block, "_refine_ls_shell.", &shell_tags, shell_rows)
}

fn software_classification_text(
    classification: cosmolkit_bio::BioSoftwareClassification,
) -> &'static str {
    // Gemmi❗✔️: std::string software_classification_to_string(SoftwareItem::Classification c) {
    // Gemmi❗✔️:   switch (c) {
    // Gemmi❗✔️:     case SoftwareItem::DataCollection: return "data collection";
    // Gemmi❗✔️:     case SoftwareItem::DataExtraction: return "data extraction";
    // Gemmi❗✔️:     case SoftwareItem::DataProcessing: return "data processing";
    // Gemmi❗✔️:     case SoftwareItem::DataReduction:  return "data reduction";
    // Gemmi❗✔️:     case SoftwareItem::DataScaling:    return "data scaling";
    // Gemmi❗✔️:     case SoftwareItem::ModelBuilding:  return "model building";
    // Gemmi❗✔️:     case SoftwareItem::Phasing:        return "phasing";
    // Gemmi❗✔️:     case SoftwareItem::Refinement:     return "refinement";
    // Gemmi❗✔️:     case SoftwareItem::Unspecified:    return "";
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   unreachable();
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    use cosmolkit_bio::BioSoftwareClassification;
    match classification {
        BioSoftwareClassification::DataCollection => "data collection",
        BioSoftwareClassification::DataExtraction => "data extraction",
        BioSoftwareClassification::DataProcessing => "data processing",
        BioSoftwareClassification::DataReduction => "data reduction",
        BioSoftwareClassification::DataScaling => "data scaling",
        BioSoftwareClassification::ModelBuilding => "model building",
        BioSoftwareClassification::Phasing => "phasing",
        BioSoftwareClassification::Refinement => "refinement",
        BioSoftwareClassification::Unspecified => "",
    }
}

fn write_tls_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.tls && st.meta.get_tls_groups() != nullptr) {
    // Gemmi❗❌:     // pdbx_refine_id doesn't make sense here, but it's required
    // Gemmi❗❌:     // by the mmCIF spec. In joint refinement, TLS constraints can't be
    // Gemmi❗❌:     // specific to a dataset, because they constrain the shared model.
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_pdbx_refine_tls.", {
    // Gemmi❗❌:         "id", "pdbx_refine_id",
    // Gemmi❗❌:         "origin_x", "origin_y", "origin_z",
    // Gemmi❗❌:         "T[1][1]", "T[2][2]", "T[3][3]", "T[1][2]", "T[1][3]", "T[2][3]",
    // Gemmi❗❌:         "L[1][1]", "L[2][2]", "L[3][3]", "L[1][2]", "L[1][3]", "L[2][3]",
    // Gemmi❗❌:         "S[1][1]", "S[1][2]", "S[1][3]",
    // Gemmi❗❌:         "S[2][1]", "S[2][2]", "S[2][3]",
    // Gemmi❗❌:         "S[3][1]", "S[3][2]", "S[3][3]"});
    // Gemmi❗❌:     for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi❗❌:       for (const TlsGroup& tls : ref.tls_groups) {
    // Gemmi❗❌:         const SMat33<double>& T = tls.T;
    // Gemmi❗❌:         const SMat33<double>& L = tls.L;
    // Gemmi❗❌:         const Mat33& S = tls.S;
    // Gemmi❗❌:         auto q = number_or_qmark;
    // Gemmi❗❌:         loop.add_row({string_or_dot(tls.id), cif::quote(ref.id),
    // Gemmi❗❌:                       q(tls.origin.x), q(tls.origin.y), q(tls.origin.z),
    // Gemmi❗❌:                       q(T.u11), q(T.u22), q(T.u33), q(T.u12), q(T.u13), q(T.u23),
    // Gemmi❗❌:                       q(L.u11), q(L.u22), q(L.u33), q(L.u12), q(L.u13), q(L.u23),
    // Gemmi❗❌:                       q(S[0][0]), q(S[0][1]), q(S[0][2]),
    // Gemmi❗❌:                       q(S[1][0]), q(S[1][1]), q(S[1][2]),
    // Gemmi❗❌:                       q(S[2][0]), q(S[2][1]), q(S[2][2])});
    // Gemmi❗❌:       }
    // Gemmi❗❌:     cif::Loop& group_loop = block.init_mmcif_loop("_pdbx_refine_tls_group.", {
    // Gemmi❗❌:         "id", "refine_tls_id", "pdbx_refine_id",
    // Gemmi❗❌:         "beg_auth_asym_id", "beg_auth_seq_id", "beg_PDB_ins_code",
    // Gemmi❗❌:         "end_auth_asym_id", "end_auth_seq_id", "end_PDB_ins_code",
    // Gemmi❗❌:         "selection_details"});
    // Gemmi❗❌:     int counter = 1;
    // Gemmi❗❌:     for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi❗❌:       for (const TlsGroup& tls : ref.tls_groups)
    // Gemmi❗❌:         for (const TlsGroup::Selection& sel : tls.selections)
    // Gemmi❗❌:           group_loop.add_row({std::to_string(counter++),
    // Gemmi❗❌:                               string_or_dot(tls.id),
    // Gemmi❗❌:                               cif::quote(ref.id),
    // Gemmi❗❌:                               string_or_qmark(sel.chain),
    // Gemmi❗❌:                               sel.res_begin.num.str(),
    // Gemmi❗❌:                               pdbx_icode(sel.res_begin),
    // Gemmi❗❌:                               string_or_qmark(sel.chain),
    // Gemmi❗❌:                               sel.res_end.num.str(),
    // Gemmi❗❌:                               pdbx_icode(sel.res_end),
    // Gemmi❗❌:                               string_or_qmark(sel.details)});
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if !groups.tls
        || !structure
            .metadata()
            .refinement
            .iter()
            .any(|refinement| !refinement.tls_groups.is_empty())
    {
        return Ok(());
    }

    let mut tls_rows = Vec::new();
    for refinement in &structure.metadata().refinement {
        for tls in &refinement.tls_groups {
            tls_rows.push(vec![
                string_or_dot(&tls.id),
                quote_cif_value(&refinement.id),
                number_or_qmark(Some(tls.origin[0])),
                number_or_qmark(Some(tls.origin[1])),
                number_or_qmark(Some(tls.origin[2])),
                number_or_qmark(Some(tls.t[0])),
                number_or_qmark(Some(tls.t[1])),
                number_or_qmark(Some(tls.t[2])),
                number_or_qmark(Some(tls.t[3])),
                number_or_qmark(Some(tls.t[4])),
                number_or_qmark(Some(tls.t[5])),
                number_or_qmark(Some(tls.l[0])),
                number_or_qmark(Some(tls.l[1])),
                number_or_qmark(Some(tls.l[2])),
                number_or_qmark(Some(tls.l[3])),
                number_or_qmark(Some(tls.l[4])),
                number_or_qmark(Some(tls.l[5])),
                number_or_qmark(Some(tls.s[0][0])),
                number_or_qmark(Some(tls.s[0][1])),
                number_or_qmark(Some(tls.s[0][2])),
                number_or_qmark(Some(tls.s[1][0])),
                number_or_qmark(Some(tls.s[1][1])),
                number_or_qmark(Some(tls.s[1][2])),
                number_or_qmark(Some(tls.s[2][0])),
                number_or_qmark(Some(tls.s[2][1])),
                number_or_qmark(Some(tls.s[2][2])),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_refine_tls.",
        &[
            "id",
            "pdbx_refine_id",
            "origin_x",
            "origin_y",
            "origin_z",
            "T[1][1]",
            "T[2][2]",
            "T[3][3]",
            "T[1][2]",
            "T[1][3]",
            "T[2][3]",
            "L[1][1]",
            "L[2][2]",
            "L[3][3]",
            "L[1][2]",
            "L[1][3]",
            "L[2][3]",
            "S[1][1]",
            "S[1][2]",
            "S[1][3]",
            "S[2][1]",
            "S[2][2]",
            "S[2][3]",
            "S[3][1]",
            "S[3][2]",
            "S[3][3]",
        ],
        tls_rows,
    )?;

    let mut selection_rows = Vec::new();
    for refinement in &structure.metadata().refinement {
        for tls in &refinement.tls_groups {
            for selection in &tls.selections {
                let begin = selection.res_begin;
                let end = selection.res_end;
                selection_rows.push(vec![
                    (selection_rows.len() + 1).to_string(),
                    string_or_dot(&tls.id),
                    quote_cif_value(&refinement.id),
                    string_or_qmark(selection.chain.as_str()),
                    seq_number_or_qmark(Some(begin)),
                    pdbx_icode_value(Some(begin)),
                    string_or_qmark(selection.chain.as_str()),
                    seq_number_or_qmark(Some(end)),
                    pdbx_icode_value(Some(end)),
                    string_or_qmark(&selection.details),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_refine_tls_group.",
        &[
            "id",
            "refine_tls_id",
            "pdbx_refine_id",
            "beg_auth_asym_id",
            "beg_auth_seq_id",
            "beg_PDB_ins_code",
            "end_auth_asym_id",
            "end_auth_seq_id",
            "end_PDB_ins_code",
            "selection_details",
        ],
        selection_rows,
    )
}

fn write_software_category(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.software && !st.meta.software.empty()) {
    // Gemmi❗❌:     bool write_all_fields = false;
    // Gemmi❗❌:     for (const SoftwareItem& item : st.meta.software)
    // Gemmi❗❌:       if (!item.date.empty() || !item.description.empty() ||
    // Gemmi❗❌:           !item.contact_author.empty() || !item.contact_author_email.empty())
    // Gemmi❗❌:         write_all_fields = true;
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_software.",
    // Gemmi❗❌:         {"pdbx_ordinal", "classification", "name", "version"});
    // Gemmi❗❌:     if (write_all_fields)
    // Gemmi❗❌:       loop.tags.insert(loop.tags.end(),
    // Gemmi❗❌:                        {"_software.date", "_software.description",
    // Gemmi❗❌:                         "_software.contact_author", "_software.contact_author_email"});
    // Gemmi❗❌:     int ordinal = 0;
    // Gemmi❗❌:     for (const SoftwareItem& item : st.meta.software) {
    // Gemmi❗❌:       loop.add_values({
    // Gemmi❗❌:           std::to_string(++ordinal),
    // Gemmi❗❌:           cif::quote(software_classification_to_string(item.classification)),
    // Gemmi❗❌:           cif::quote(item.name),
    // Gemmi❗❌:           string_or_dot(item.version)});
    // Gemmi❗❌:       if (write_all_fields)
    // Gemmi❗❌:         loop.add_values({
    // Gemmi❗❌:             string_or_qmark(item.date),
    // Gemmi❗❌:             string_or_qmark(item.description),
    // Gemmi❗❌:             string_or_qmark(item.contact_author),
    // Gemmi❗❌:             string_or_qmark(item.contact_author_email)});
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if !groups.software || structure.metadata().software.is_empty() {
        return Ok(());
    }
    let write_all_fields = structure.metadata().software.iter().any(|item| {
        !item.date.is_empty()
            || !item.description.is_empty()
            || !item.contact_author.is_empty()
            || !item.contact_author_email.is_empty()
    });
    let mut tags = vec!["pdbx_ordinal", "classification", "name", "version"];
    if write_all_fields {
        tags.extend([
            "date",
            "description",
            "contact_author",
            "contact_author_email",
        ]);
    }
    let mut rows = Vec::with_capacity(structure.metadata().software.len());
    for (index, item) in structure.metadata().software.iter().enumerate() {
        let mut row = vec![
            (index + 1).to_string(),
            quote_cif_value(software_classification_text(item.classification)),
            quote_cif_value(&item.name),
            string_or_dot(&item.version),
        ];
        if write_all_fields {
            row.extend([
                string_or_qmark(&item.date),
                string_or_qmark(&item.description),
                string_or_qmark(&item.contact_author),
                string_or_qmark(&item.contact_author_email),
            ]);
        }
        rows.push(row);
    }
    add_mmcif_rows(block, "_software.", &tags, rows)
}

fn add_mmcif_rows(
    block: &mut CifBlock,
    category: &str,
    suffixes: &[&str],
    rows: Vec<Vec<String>>,
) -> Result<(), BioMmcifWriteError> {
    // Behavior: the canonical CIF owner checks each row and owns category replacement.
    // Cost: buffers rows before appending; extra row-vector allocations versus source direct append.
    let loop_ = block.init_mmcif_loop(category, suffixes)?;
    for row in rows {
        loop_.add_row(row)?;
    }
    Ok(())
}
fn get_number_obs(refinement: &BioBasicRefinementInfo) -> Option<i32> {
    // Gemmi❗✔️: int get_number_obs(const BasicRefinementInfo& ref) {
    // Gemmi❗✔️:   int nobs = ref.reflection_count;
    // Gemmi❗✔️:   if (nobs == -1 && ref.rfree_set_count >= 0 && ref.work_set_count >= 0)
    // Gemmi❗✔️:     nobs = ref.work_set_count + ref.rfree_set_count;
    // Gemmi❗✔️:   return nobs;
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    let mut count = refinement.reflection_count;
    if count == -1 && refinement.rfree_set_count >= 0 && refinement.work_set_count >= 0 {
        count = refinement
            .work_set_count
            .wrapping_add(refinement.rfree_set_count);
    }
    Some(count)
}
fn get_number_work(refinement: &BioBasicRefinementInfo) -> Option<i32> {
    // Gemmi❗✔️: int get_number_work(const BasicRefinementInfo& ref) {
    // Gemmi❗✔️:   int nwork = ref.work_set_count;
    // Gemmi❗✔️:   if (nwork == -1 && ref.rfree_set_count >= 0 && ref.reflection_count >= 0)
    // Gemmi❗✔️:     nwork = ref.reflection_count - ref.rfree_set_count;
    // Gemmi❗✔️:   return nwork;
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    let mut count = refinement.work_set_count;
    if count == -1 && refinement.rfree_set_count >= 0 && refinement.reflection_count >= 0 {
        count = refinement
            .reflection_count
            .wrapping_sub(refinement.rfree_set_count);
    }
    Some(count)
}
fn pdbx_icode_value(value: Option<PdbSeqId>) -> String {
    // Gemmi❗✔️: inline std::string pdbx_icode(const SeqId& seqid) {
    // Gemmi❗✔️:   return std::string(1, seqid.has_icode() ? seqid.icode : '?');
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    value
        .and_then(|seq| seq.ins_code())
        .map_or_else(|| "?".to_owned(), |b| char::from(b).to_string())
}

fn quote_cif_value(value: impl AsRef<str>) -> String {
    crate::cif::quote_cif_value(value.as_ref().to_owned())
}

pub(super) fn write_secondary_mmcif_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
    entry_id: &str,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.ncs)
    // Gemmi❗❌:     write_ncs_oper(st, block);
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.struct_asym) {
    // Gemmi❗❌:     cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
    // Gemmi❗❌:                                                  {"id", "entity_id"});
    // Gemmi❗❌:     for (const Chain& chain : st.models[0].chains)
    // Gemmi❗❌:       for (ConstResidueSpan& sub : chain.subchains()) {
    // Gemmi❗❌:         const std::string& sub_id = sub.subchain_id();
    // Gemmi❗❌:         if (!sub_id.empty()) {
    // Gemmi❗❌:           const Entity* ent = find_entity_of_subchain(sub_id, st.entities);
    // Gemmi❗❌:           asym_loop.add_row({sub_id, (ent ? qchain(ent->name) : "?")});
    // Gemmi❗❌:         }
    // Gemmi❗❌:       }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   bool nontrivial_origx = st.has_origx && !st.origx.is_identity();
    // Gemmi❗❌:   if (groups.origx && nontrivial_origx) { // _database_PDB_matrix (ORIGX)
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_database_PDB_matrix.");
    // Gemmi❗❌:     span.set_pair("_database_PDB_matrix.entry_id", id);
    // Gemmi❗❌:     std::string tag_mat = "_database_PDB_matrix.origx[0][0]";
    // Gemmi❗❌:     std::string tag_vec = "_database_PDB_matrix.origx_vector[0]";
    // Gemmi❗❌:     for (int i = 0; i < 3; ++i) {
    // Gemmi❗❌:       tag_mat[27] += 1;  // origx[0] -> origx[1] -> origx[2]
    // Gemmi❗❌:       tag_vec[34] += 1;
    // Gemmi❗❌:       for (int j = 0; j < 3; ++j) {
    // Gemmi❗❌:         tag_mat[30] = '1' + j;
    // Gemmi❗❌:         span.set_pair(tag_mat, to_str(st.origx.mat[i][j]));
    // Gemmi❗❌:       }
    // Gemmi❗❌:       span.set_pair(tag_vec, to_str(st.origx.vec.at(i)));
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.struct_conf && !st.helices.empty()) {
    // Gemmi❗❌:     cif::Loop& struct_conf_loop = block.init_mmcif_loop("_struct_conf.",
    // Gemmi❗❌:         {"conf_type_id", "id",
    // Gemmi❗❌:          "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
    // Gemmi❗❌:          "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
    // Gemmi❗❌:          "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
    // Gemmi❗❌:          "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code",
    // Gemmi❗❌:          "pdbx_PDB_helix_class", "pdbx_PDB_helix_length"});
    // Gemmi❗❌:     int count = 0;
    // Gemmi❗❌:     for (const Helix& helix : st.helices) {
    // Gemmi❗❌:       const_CRA cra1 = st.models[0].find_cra(helix.start);
    // Gemmi❗❌:       const_CRA cra2 = st.models[0].find_cra(helix.end);
    // Gemmi❗❌:       if (!cra1.residue || !cra2.residue)
    // Gemmi❗❌:         continue;
    // Gemmi❗❌:       struct_conf_loop.add_row({
    // Gemmi❗❌:         "HELX_P",                                    // conf_type_id
    // Gemmi❗❌:         "H" + std::to_string(++count),               // id
    // Gemmi❗❌:         qchain(cra1.chain->name),                    // beg_auth_asym_id
    // Gemmi❗❌:         subchain_or_dot(*cra1.residue),              // beg_label_asym_id
    // Gemmi❗❌:         cra1.residue->name,                          // beg_label_comp_id
    // Gemmi❗❌:         cra1.residue->label_seq.str(),               // beg_label_seq_id
    // Gemmi❗❌:         cra1.residue->seqid.num.str(),               // beg_auth_seq_id
    // Gemmi❗❌:         pdbx_icode(*cra1.residue),                   // beg_PDB_ins_code
    // Gemmi❗❌:         qchain(cra2.chain->name),                    // end_auth_asym_id
    // Gemmi❗❌:         subchain_or_dot(*cra2.residue),              // end_label_asym_id
    // Gemmi❗❌:         cra2.residue->name,                          // end_label_comp_id
    // Gemmi❗❌:         cra2.residue->label_seq.str(),               // end_label_seq_id
    // Gemmi❗❌:         cra2.residue->seqid.num.str(),               // end_auth_seq_id
    // Gemmi❗❌:         pdbx_icode(*cra2.residue),                   // end_PDB_ins_code
    // Gemmi❗❌:         std::to_string((int)helix.pdb_helix_class),  // pdbx_PDB_helix_class
    // Gemmi❗❌:         int_or_qmark(helix.length)                   // pdbx_PDB_helix_length
    // Gemmi❗❌:       });
    // Gemmi❗❌:     }
    // Gemmi❗❌:     if (count != 0)
    // Gemmi❗❌:       block.set_pair("_struct_conf_type.id", "HELX_P");
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if structure.models().is_empty() {
        return Ok(());
    }

    if groups.ncs {
        super::ncs::write_ncs_oper(structure, block)?;
    }
    if groups.struct_asym {
        let mut rows = Vec::new();
        for chain in structure.models()[0]
            .chain_span()
            .slice(structure.chains())?
        {
            let mut previous = None;
            for residue in chain.residue_span().slice(structure.residues())? {
                let subchain = residue.source().subchain_id().unwrap_or("");
                if previous == Some(subchain) {
                    continue;
                }
                previous = Some(subchain);
                if subchain.is_empty() {
                    continue;
                }
                let entity_id = structure.find_entity_of_subchain(subchain).map_or_else(
                    || "?".to_owned(),
                    |(_, entity)| quote_cif_value(entity.source().source_entity_id()),
                );
                rows.push(vec![subchain.to_owned(), entity_id]);
            }
        }
        add_mmcif_rows(block, "_struct_asym.", &["id", "entity_id"], rows)?;
    }

    let nontrivial_origx = structure.source_state().has_origx
        && !transform_is_exact_identity(&structure.source_state().origx);
    if groups.origx {
        super::origx::write_origx_category(
            structure.source_state().has_origx,
            &structure.source_state().origx,
            entry_id,
            block,
        );
    }
    if groups.struct_conf && !structure.helices().is_empty() {
        let mut rows = Vec::new();
        for helix in structure.helices() {
            let Some(begin) = find_cra(structure, 0, &helix.start)? else {
                continue;
            };
            let Some(end) = find_cra(structure, 0, &helix.end)? else {
                continue;
            };
            let begin_values = secondary_residue_values(structure, begin, &helix.start)?;
            let end_values = secondary_residue_values(structure, end, &helix.end)?;
            let count = rows.len() + 1;
            rows.push(vec![
                "HELX_P".to_string(),
                format!("H{count}"),
                begin_values[0].clone(),
                begin_values[1].clone(),
                begin_values[2].clone(),
                begin_values[3].clone(),
                begin_values[4].clone(),
                begin_values[5].clone(),
                end_values[0].clone(),
                end_values[1].clone(),
                end_values[2].clone(),
                end_values[3].clone(),
                end_values[4].clone(),
                end_values[5].clone(),
                (helix.pdb_helix_class as i32).to_string(),
                if helix.length == -1 {
                    "?".to_string()
                } else {
                    helix.length.to_string()
                },
            ]);
        }
        let has_rows = !rows.is_empty();
        add_mmcif_rows(
            block,
            "_struct_conf.",
            &[
                "conf_type_id",
                "id",
                "beg_auth_asym_id",
                "beg_label_asym_id",
                "beg_label_comp_id",
                "beg_label_seq_id",
                "beg_auth_seq_id",
                "pdbx_beg_PDB_ins_code",
                "end_auth_asym_id",
                "end_label_asym_id",
                "end_label_comp_id",
                "end_label_seq_id",
                "end_auth_seq_id",
                "pdbx_end_PDB_ins_code",
                "pdbx_PDB_helix_class",
                "pdbx_PDB_helix_length",
            ],
            rows,
        )?;
        if has_rows {
            block.set_pair_in_category(None, "_struct_conf_type.id", "HELX_P".to_owned());
        }
    }

    write_secondary_sheet_categories(structure, block, groups)?;
    write_secondary_tail_categories(structure, block, groups, entry_id, nontrivial_origx)
}

fn write_secondary_sheet_categories(
    structure: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   // _struct_sheet*
    // Gemmi❗❌:   if (groups.struct_sheet && !st.sheets.empty()) {
    // Gemmi❗❌:     cif::Loop& sheet_loop = block.init_mmcif_loop("_struct_sheet.",
    // Gemmi❗❌:                                                   {"id", "number_strands"});
    // Gemmi❗❌:     for (const Sheet& sheet : st.sheets)
    // Gemmi❗❌:       sheet_loop.add_row({string_or_dot(sheet.name),
    // Gemmi❗❌:                           std::to_string(sheet.strands.size())});
    // Gemmi❗❌:
    // Gemmi❗❌:     cif::Loop& order_loop = block.init_mmcif_loop("_struct_sheet_order.",
    // Gemmi❗❌:                     {"sheet_id", "range_id_1", "range_id_2", "sense"});
    // Gemmi❗❌:     for (const Sheet& sheet : st.sheets)
    // Gemmi❗❌:       for (size_t i = 1; i < sheet.strands.size(); ++i) {
    // Gemmi❗❌:         const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi❗❌:         if (strand.sense != 0)
    // Gemmi❗❌:           order_loop.add_row({string_or_dot(sheet.name),
    // Gemmi❗❌:                               std::to_string(i), std::to_string(i+1),
    // Gemmi❗❌:                               strand.sense > 0 ? "parallel" : "anti-parallel"});
    // Gemmi❗❌:       }
    // Gemmi❗❌:
    // Gemmi❗❌:     cif::Loop& range_loop = block.init_mmcif_loop("_struct_sheet_range.",
    // Gemmi❗❌:         {"sheet_id", "id",
    // Gemmi❗❌:          "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
    // Gemmi❗❌:          "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
    // Gemmi❗❌:          "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
    // Gemmi❗❌:          "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code"});
    // Gemmi❗❌:     for (const Sheet& sheet : st.sheets)
    // Gemmi❗❌:       for (size_t i = 0; i < sheet.strands.size(); ++i) {
    // Gemmi❗❌:         const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi❗❌:         const_CRA cra1 = st.models[0].find_cra(strand.start);
    // Gemmi❗❌:         const_CRA cra2 = st.models[0].find_cra(strand.end);
    // Gemmi❗❌:         if (!cra1.residue || !cra2.residue)
    // Gemmi❗❌:           continue;
    // Gemmi❗❌:         range_loop.add_row({
    // Gemmi❗❌:           string_or_dot(sheet.name),            // sheet_id
    // Gemmi❗❌:           std::to_string(i+1),                  // id
    // Gemmi❗❌:           qchain(cra1.chain->name),             // beg_auth_asym_id
    // Gemmi❗❌:           subchain_or_dot(*cra1.residue),       // beg_label_asym_id
    // Gemmi❗❌:           cra1.residue->name,                   // beg_label_comp_id
    // Gemmi❗❌:           cra1.residue->label_seq.str(),        // beg_label_seq_id
    // Gemmi❗❌:           cra1.residue->seqid.num.str(),        // beg_auth_seq_id
    // Gemmi❗❌:           pdbx_icode(*cra1.residue),            // beg_PDB_ins_code
    // Gemmi❗❌:           qchain(cra2.chain->name),             // end_auth_asym_id
    // Gemmi❗❌:           subchain_or_dot(*cra2.residue),       // end_label_asym_id
    // Gemmi❗❌:           cra2.residue->name,                   // end_label_comp_id
    // Gemmi❗❌:           cra2.residue->label_seq.str(),        // end_label_seq_id
    // Gemmi❗❌:           cra2.residue->seqid.num.str(),        // end_auth_seq_id
    // Gemmi❗❌:           pdbx_icode(*cra2.residue)             // end_PDB_ins_code
    // Gemmi❗❌:         });
    // Gemmi❗❌:     }
    // Gemmi❗❌:
    // Gemmi❗❌:     cif::Loop& hbond_loop = block.init_mmcif_loop("_pdbx_struct_sheet_hbond.",
    // Gemmi❗❌:         {"sheet_id", "range_id_1", "range_id_2",
    // Gemmi❗❌:          "range_1_auth_asym_id", "range_1_label_asym_id",
    // Gemmi❗❌:          "range_1_label_comp_id", "range_1_label_seq_id", "range_1_auth_seq_id",
    // Gemmi❗❌:          "range_1_PDB_ins_code", "range_1_label_atom_id",
    // Gemmi❗❌:          "range_2_auth_asym_id", "range_2_label_asym_id",
    // Gemmi❗❌:          "range_2_label_comp_id", "range_2_label_seq_id", "range_2_auth_seq_id",
    // Gemmi❗❌:          "range_2_PDB_ins_code", "range_2_label_atom_id"});
    // Gemmi❗❌:     for (const Sheet& sheet : st.sheets)
    // Gemmi❗❌:       for (size_t i = 1; i < sheet.strands.size(); ++i) {
    // Gemmi❗❌:         const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi❗❌:         if (strand.hbond_atom2.atom_name.empty())
    // Gemmi❗❌:           continue;
    // Gemmi❗❌:         // hbond_atomN is not a full atom "address": altloc is missing
    // Gemmi❗❌:         const_CRA cra1 = st.models[0].find_cra(strand.hbond_atom1);
    // Gemmi❗❌:         const_CRA cra2 = st.models[0].find_cra(strand.hbond_atom2);
    // Gemmi❗❌:         if (!cra1.residue || !cra2.residue)
    // Gemmi❗❌:           continue;
    // Gemmi❗❌:         hbond_loop.add_row({
    // Gemmi❗❌:           string_or_dot(sheet.name),                  // sheet_id
    // Gemmi❗❌:           std::to_string(i),                          // range_id_1
    // Gemmi❗❌:           std::to_string(i+1),                        // range_id_2
    // Gemmi❗❌:           qchain(cra1.chain->name),                   // range_1_auth_asym_id
    // Gemmi❗❌:           subchain_or_dot(*cra1.residue),             // range_1_label_asym_id
    // Gemmi❗❌:           cra1.residue->name,                         // range_1_label_comp_id
    // Gemmi❗❌:           cra1.residue->label_seq.str(),              // range_1_label_seq_id
    // Gemmi❗❌:           cra1.residue->seqid.num.str(),              // range_1_auth_seq_id
    // Gemmi❗❌:           pdbx_icode(*cra1.residue),                  // range_1_PDB_ins_code
    // Gemmi❗❌:           cif::quote(strand.hbond_atom1.atom_name),   // range_1_label_atom_id
    // Gemmi❗❌:           qchain(cra2.chain->name),                   // range_2_auth_asym_id
    // Gemmi❗❌:           subchain_or_dot(*cra2.residue),             // range_2_label_asym_id
    // Gemmi❗❌:           cra2.residue->name,                         // range_2_label_comp_id
    // Gemmi❗❌:           cra2.residue->label_seq.str(),              // range_2_label_seq_id
    // Gemmi❗❌:           cra2.residue->seqid.num.str(),              // range_2_auth_seq_id
    // Gemmi❗❌:           pdbx_icode(*cra2.residue),                  // range_2_PDB_ins_code
    // Gemmi❗❌:           cif::quote(strand.hbond_atom2.atom_name)    // range_2_label_atom_id
    // Gemmi❗❌:         });
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    if !groups.struct_sheet || structure.sheets().is_empty() {
        return Ok(());
    }

    add_mmcif_rows(
        block,
        "_struct_sheet.",
        &["id", "number_strands"],
        structure
            .sheets()
            .iter()
            .map(|sheet| vec![string_or_dot(&sheet.name), sheet.strands.len().to_string()])
            .collect(),
    )?;

    let mut order_rows = Vec::new();
    for sheet in structure.sheets() {
        for (index, strand) in sheet.strands.iter().enumerate().skip(1) {
            if strand.sense != 0 {
                order_rows.push(vec![
                    string_or_dot(&sheet.name),
                    index.to_string(),
                    (index + 1).to_string(),
                    if strand.sense > 0 {
                        "parallel"
                    } else {
                        "anti-parallel"
                    }
                    .to_string(),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_struct_sheet_order.",
        &["sheet_id", "range_id_1", "range_id_2", "sense"],
        order_rows,
    )?;

    let mut range_rows = Vec::new();
    for sheet in structure.sheets() {
        for (index, strand) in sheet.strands.iter().enumerate() {
            let Some(begin) = find_cra(structure, 0, &strand.start)? else {
                continue;
            };
            let Some(end) = find_cra(structure, 0, &strand.end)? else {
                continue;
            };
            let begin_values = secondary_residue_values(structure, begin, &strand.start)?;
            let end_values = secondary_residue_values(structure, end, &strand.end)?;
            range_rows.push(vec![
                string_or_dot(&sheet.name),
                (index + 1).to_string(),
                begin_values[0].clone(),
                begin_values[1].clone(),
                begin_values[2].clone(),
                begin_values[3].clone(),
                begin_values[4].clone(),
                begin_values[5].clone(),
                end_values[0].clone(),
                end_values[1].clone(),
                end_values[2].clone(),
                end_values[3].clone(),
                end_values[4].clone(),
                end_values[5].clone(),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_struct_sheet_range.",
        &[
            "sheet_id",
            "id",
            "beg_auth_asym_id",
            "beg_label_asym_id",
            "beg_label_comp_id",
            "beg_label_seq_id",
            "beg_auth_seq_id",
            "pdbx_beg_PDB_ins_code",
            "end_auth_asym_id",
            "end_label_asym_id",
            "end_label_comp_id",
            "end_label_seq_id",
            "end_auth_seq_id",
            "pdbx_end_PDB_ins_code",
        ],
        range_rows,
    )?;

    let mut hbond_rows = Vec::new();
    for sheet in structure.sheets() {
        for (index, strand) in sheet.strands.iter().enumerate().skip(1) {
            if strand.hbond_atom2.logical_atom_name().is_empty() {
                continue;
            }
            let Some(left) = find_cra(structure, 0, &strand.hbond_atom1)? else {
                continue;
            };
            let Some(right) = find_cra(structure, 0, &strand.hbond_atom2)? else {
                continue;
            };
            let left_values = secondary_residue_values(structure, left, &strand.hbond_atom1)?;
            let right_values = secondary_residue_values(structure, right, &strand.hbond_atom2)?;
            hbond_rows.push(vec![
                string_or_dot(&sheet.name),
                index.to_string(),
                (index + 1).to_string(),
                left_values[0].clone(),
                left_values[1].clone(),
                left_values[2].clone(),
                left_values[3].clone(),
                left_values[4].clone(),
                left_values[5].clone(),
                quote_cif_value(&strand.hbond_atom1.logical_atom_name()),
                right_values[0].clone(),
                right_values[1].clone(),
                right_values[2].clone(),
                right_values[3].clone(),
                right_values[4].clone(),
                right_values[5].clone(),
                quote_cif_value(&strand.hbond_atom2.logical_atom_name()),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_struct_sheet_hbond.",
        &[
            "sheet_id",
            "range_id_1",
            "range_id_2",
            "range_1_auth_asym_id",
            "range_1_label_asym_id",
            "range_1_label_comp_id",
            "range_1_label_seq_id",
            "range_1_auth_seq_id",
            "range_1_PDB_ins_code",
            "range_1_label_atom_id",
            "range_2_auth_asym_id",
            "range_2_label_asym_id",
            "range_2_label_comp_id",
            "range_2_label_seq_id",
            "range_2_auth_seq_id",
            "range_2_PDB_ins_code",
            "range_2_label_atom_id",
        ],
        hbond_rows,
    )?;
    Ok(())
}

#[derive(Clone, Copy)]
struct WriterCra {
    residue_index: usize,
    atom_index: Option<usize>,
}
fn find_cra(
    data: &BioStructureData,
    model: usize,
    address: &AtomAddress,
) -> Result<Option<WriterCra>, BioMmcifWriteError> {
    // Gemmi❗✔️:     const_CRA cra1 = st.models[0].find_cra(helix.start);
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    // IO keeps only output indices; BIO owns source address matching.
    Ok(data
        .find_cra(BioModelId::new(model as u32), address, false)?
        .map(|(_, residue, atom)| WriterCra {
            residue_index: residue.index(),
            atom_index: atom.map(|id| id.index()),
        }))
}
fn secondary_residue_values(
    data: &BioStructureData,
    cra: WriterCra,
    _address: &AtomAddress,
) -> Result<[String; 6], BioMmcifWriteError> {
    // Gemmi❗❌:         qchain(cra1.chain->name),                    // beg_auth_asym_id
    // Gemmi❗❌:         subchain_or_dot(*cra1.residue),              // beg_label_asym_id
    // Gemmi❗❌:         cra1.residue->name,                          // beg_label_comp_id
    // Gemmi❗❌:         cra1.residue->label_seq.str(),               // beg_label_seq_id
    // Gemmi❗❌:         cra1.residue->seqid.num.str(),               // beg_auth_seq_id
    // Gemmi❗❌:         pdbx_icode(*cra1.residue),                   // beg_PDB_ins_code
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let row = &data.residues()[cra.residue_index];
    let chain = &data.chains()[row.chain_id().index()];
    Ok([
        quote_cif_value(
            chain
                .source()
                .auth_chain_id()
                .map_or_else(String::new, |id| id.as_str().to_owned()),
        ),
        string_or_dot(row.source().subchain_id().unwrap_or("")),
        residue_name_logical_view(&row.name(), data.input_format()).to_owned(),
        row.source()
            .label_seq_id()
            .map_or_else(|| "?".to_owned(), |n| n.to_string()),
        seq_number_or_qmark(row.source().seq_id()),
        pdbx_icode_value(row.source().seq_id()),
    ])
}
fn seq_number_or_qmark(seq: Option<PdbSeqId>) -> String {
    seq.filter(|s| s.seq_num() != i32::MIN)
        .map_or_else(|| "?".to_owned(), |s| s.seq_num().to_string())
}
fn transform_is_exact_identity(t: &BioTransform) -> bool {
    !super::origx::origx_is_nontrivial(true, t)
}
fn transform_approx(a: &BioTransform, b: &BioTransform) -> bool {
    a.approx(b, 1e-9)
}
fn first_model_subchain_to_chain(
    data: &BioStructureData,
) -> Result<std::collections::BTreeMap<&str, String>, BioMmcifWriteError> {
    // Gemmi❗❌:   std::map<std::string, std::string> subchain_to_chain() const {
    // Gemmi❗❌:     std::map<std::string, std::string> mapping;
    // Gemmi❗❌:     for (const Chain& chain : chains) {
    // Gemmi❗❌:       std::string prev;
    // Gemmi❗❌:       for (const Residue& res : chain.residues)
    // Gemmi❗❌:         if (!res.subchain.empty() && res.subchain != prev) {
    // Gemmi❗❌:           prev = res.subchain;
    // Gemmi❗❌:           mapping[res.subchain] = chain.name;
    // Gemmi❗❌:         }
    // Gemmi❗❌:     }
    // Gemmi❗❌:     return mapping;
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let mut map = std::collections::BTreeMap::new();
    if let Some(model) = data.models().first() {
        for chain in model.chain_span().slice(data.chains())? {
            let chain_name = chain
                .source()
                .auth_chain_id()
                .map_or_else(String::new, |id| id.as_str().to_owned());
            let mut previous = "";
            for residue in chain.residue_span().slice(data.residues())? {
                let sub = residue.source().subchain_id().unwrap_or("");
                if !sub.is_empty() && sub != previous {
                    previous = sub;
                    map.insert(sub, chain_name.clone());
                }
            }
        }
    }
    Ok(map)
}
pub(super) fn write_primary_mmcif_categories(
    data: &BioStructureData,
    block: &mut CifBlock,
    groups: &BioMmcifWriteParams,
) -> Result<String, BioMmcifWriteError> {
    // Gemmi❗❌: void update_mmcif_block(const Structure& st, cif::Block& block, MmcifOutputGroups groups) {
    // Gemmi❗❌:   if (st.models.empty())
    // Gemmi❗❌:     return;
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.block_name)
    // Gemmi❗❌:     block.name = is_valid_block_name(st.name) ? st.name : "model";
    // Gemmi❗❌:
    // Gemmi❗❌:   auto e_id = st.info.find("_entry.id");
    // Gemmi❗❌:   std::string id = cif::quote(e_id != st.info.end() ? e_id->second : block.name);
    // Gemmi❗❌:   if (groups.entry)
    // Gemmi❗❌:     block.set_pair("_entry.id", id);
    // Gemmi❗❌:   else if (const std::string* val = block.find_value("_entry.id"))
    // Gemmi❗❌:     id = *val;
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.database_status) {
    // Gemmi❗❌:     auto initial_date = st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
    // Gemmi❗❌:     if (initial_date != st.info.end() && !initial_date->second.empty()) {
    // Gemmi❗❌:       cif::ItemSpan span(block.items, "_pdbx_database_status.");
    // Gemmi❗❌:       span.set_pair("_pdbx_database_status.entry_id", id);
    // Gemmi❗❌:       span.set_pair(initial_date->first, initial_date->second);
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.author && !st.meta.authors.empty()) {
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_audit_author.", {"pdbx_ordinal", "name"});
    // Gemmi❗❌:     int n = 0;
    // Gemmi❗❌:     for (const std::string& author : st.meta.authors)
    // Gemmi❗❌:       loop.add_row({std::to_string(++n), cif::quote(author)});
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.cell) {
    // Gemmi❗❌:     cif::ItemSpan cell_span(block.items, "_cell.");
    // Gemmi❗❌:     cell_span.set_pair("_cell.entry_id", id);
    // Gemmi❗❌:     write_cell_parameters(st.cell, cell_span);
    // Gemmi❗❌:     auto z_pdb = st.info.find("_cell.Z_PDB");
    // Gemmi❗❌:     if (z_pdb != st.info.end())
    // Gemmi❗❌:       cell_span.set_pair(z_pdb->first, z_pdb->second);
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.symmetry) {
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_symmetry.");
    // Gemmi❗❌:     span.set_pair("_symmetry.entry_id", id);
    // Gemmi❗❌:     span.set_pair("_symmetry.space_group_name_H-M",
    // Gemmi❗❌:                    cif::quote(st.spacegroup_hm));
    // Gemmi❗❌:     if (const SpaceGroup* sg = st.find_spacegroup())
    // Gemmi❗❌:       span.set_pair("_symmetry.Int_Tables_number", std::to_string(sg->number));
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.entity) {
    // Gemmi❗❌:     cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
    // Gemmi❗❌:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       entity_loop.add_row({qchain(ent.name),
    // Gemmi❗❌:                            entity_type_to_string(ent.entity_type)});
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   std::map<std::string, std::string> subs_to_strands;
    // Gemmi❗❌:   if (groups.entity_poly || groups.struct_ref)
    // Gemmi❗❌:     subs_to_strands = st.models[0].subchain_to_chain();
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.entity_poly) {
    // Gemmi❗❌:     // If the _entity_poly category is included when depositing to the PDB,
    // Gemmi❗❌:     // it must contain entity_id, type, pdbx_seq_one_letter_code
    // Gemmi❗❌:     // and pdbx_strand_id. The last one is not documented as required,
    // Gemmi❗❌:     // but OneDep shows error when it's not included.
    // Gemmi❗❌:     cif::Loop& ent_poly_loop = block.init_mmcif_loop("_entity_poly.",
    // Gemmi❗❌:         {"entity_id", "type", "pdbx_strand_id", "pdbx_seq_one_letter_code"});
    // Gemmi❗❌:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       if (ent.entity_type == EntityType::Polymer) {
    // Gemmi❗❌:         if (ent.polymer_type == PolymerType::Unknown)
    // Gemmi❗❌:           continue;  // not sure what to do here
    // Gemmi❗❌:         ResidueKind kind = sequence_kind(ent.polymer_type);
    // Gemmi❗❌:         std::string seq1 = pdbx_one_letter_code(ent.full_sequence, kind);
    // Gemmi❗❌:         std::string strand_ids;
    // Gemmi❗❌:         for (const std::string& sub : ent.subchains) {
    // Gemmi❗❌:           auto strand_id = subs_to_strands.find(sub);
    // Gemmi❗❌:           if (strand_id != subs_to_strands.end()) {
    // Gemmi❗❌:             if (!strand_ids.empty())
    // Gemmi❗❌:               strand_ids += ',';
    // Gemmi❗❌:             strand_ids += strand_id->second;
    // Gemmi❗❌:           }
    // Gemmi❗❌:         }
    // Gemmi❗❌:         ent_poly_loop.add_row({qchain(ent.name),
    // Gemmi❗❌:                                polymer_type_to_string(ent.polymer_type),
    // Gemmi❗❌:                                string_or_qmark(strand_ids),
    // Gemmi❗❌:                                string_or_qmark(seq1)});
    // Gemmi❗❌:       }
    // Gemmi❗❌:   }
    // Gemmi❗❌:   if (groups.chem_comp) {
    // Gemmi❗❌:     std::set<std::string> resnames;
    // Gemmi❗❌:     for (const Model& model : st.models)
    // Gemmi❗❌:       for (const Chain& chain : model.chains)
    // Gemmi❗❌:         for (const Residue& res : chain.residues)
    // Gemmi❗❌:           resnames.insert(res.name);
    // Gemmi❗❌:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       for (const std::string& item : ent.full_sequence)
    // Gemmi❗❌:         resnames.insert(Entity::first_mon(item));
    // Gemmi❗❌:     cif::Loop& chem_comp_loop = block.init_mmcif_loop("_chem_comp.", {"id", "type"});
    // Gemmi❌❌:     if (!st.shortened_ccd_codes.empty())
    // Gemmi❌❌:       chem_comp_loop.tags.push_back("_chem_comp.three_letter_code");
    // Gemmi❗❌:     for (const std::string& name : resnames) {
    // Gemmi❗❌:       chem_comp_loop.values.push_back(cif::quote(name));
    // Gemmi❗❌:       chem_comp_loop.values.push_back(".");
    // Gemmi❌❌:       if (!st.shortened_ccd_codes.empty()) {
    // Gemmi❌❌:         chem_comp_loop.values.push_back(cif::quote(name));
    // Gemmi❌❌:         for (const auto& old_new : st.shortened_ccd_codes)
    // Gemmi❌❌:           if (old_new.second == name)
    // Gemmi❌❌:             chem_comp_loop.values.back() = old_new.first;
    // Gemmi❌❌:       }
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let state = data.source_state();
    if groups.block_name {
        let name = &state.name;
        block.set_name(
            if !name.is_empty() && name.bytes().all(|b| (b'!'..=b'~').contains(&b)) {
                name.clone()
            } else {
                "model".to_owned()
            },
        );
    }
    let mut id = quote_cif_value(
        state
            .info
            .get("_entry.id")
            .map_or(block.name(), String::as_str),
    );
    if groups.entry {
        block.set_pair_in_category(None, "_entry.id", id.clone());
    } else if let Some(value) = block.find_value("_entry.id") {
        id = value.raw().to_owned();
    }
    if groups.database_status {
        super::database_status::write_database_status_category(&state.info, &id, block);
    }
    if groups.author {
        super::author::write_author_category(&data.metadata().authors, block)?;
    }
    if groups.cell {
        let default_crystal;
        let crystal = match data.crystal() {
            Some(crystal) => crystal,
            None => {
                default_crystal = self::default_crystal();
                &default_crystal
            }
        };
        super::cell_category::write_cell_category(crystal, &id, block);
    }
    if groups.symmetry {
        block.set_pair_in_category(Some("_symmetry."), "_symmetry.entry_id", id.clone());
        block.set_pair_in_category(
            Some("_symmetry."),
            "_symmetry.space_group_name_H-M",
            quote_cif_value(
                data.crystal()
                    .and_then(|c| c.space_group_hm())
                    .unwrap_or(""),
            ),
        );
        if let Some(number) = data.crystal().and_then(|c| c.space_group_number()) {
            block.set_pair_in_category(
                Some("_symmetry."),
                "_symmetry.Int_Tables_number",
                number.to_string(),
            );
        }
    }
    if groups.entity {
        super::entity::write_entity_category(data.entities(), block)?;
    }
    let subchain_map = if groups.entity_poly || groups.struct_ref {
        first_model_subchain_to_chain(data)?
    } else {
        std::collections::BTreeMap::new()
    };
    if groups.entity_poly {
        let mut rows = Vec::new();
        for entity in data.entities() {
            if entity.kind() != EntityKind::Polymer {
                continue;
            }
            let Some(polymer) = polymer_kind_text(entity.polymer_kind()) else {
                continue;
            };
            let strand_ids = entity
                .subchains()
                .iter()
                .filter_map(|sub| subchain_map.get(sub.as_str()).map(String::as_str))
                .collect::<Vec<_>>()
                .join(",");
            rows.push(vec![
                quote_cif_value(entity.source().source_entity_id()),
                polymer.to_owned(),
                string_or_qmark(&strand_ids),
                string_or_qmark(&pdbx_one_letter_code(
                    entity.full_sequence(),
                    entity.polymer_kind(),
                )),
            ]);
        }
        add_mmcif_rows(
            block,
            "_entity_poly.",
            &[
                "entity_id",
                "type",
                "pdbx_strand_id",
                "pdbx_seq_one_letter_code",
            ],
            rows,
        )?;
    }
    if groups.struct_ref {
        write_struct_ref_categories(data, block, &subchain_map, &id)?;
    }
    if groups.chem_comp {
        let mut names = std::collections::BTreeSet::new();
        for row in data.residues() {
            names.insert(residue_name_logical_view(&row.name(), data.input_format()).to_owned());
        }
        for entity in data.entities() {
            for item in entity.full_sequence() {
                names.insert(BioEntityRow::first_mon(item).to_owned());
            }
        }
        // No shortened-CCD mutation API exists in the modeled BIO data; Gemmi's read defaults keep its mapping empty.
        add_mmcif_rows(
            block,
            "_chem_comp.",
            &["id", "type"],
            names
                .into_iter()
                .map(|name| vec![quote_cif_value(name), ".".to_owned()])
                .collect(),
        )?;
    }
    write_experiment_categories(data, block, *groups, &id)?;
    write_reflection_categories(data, block, *groups, &id)?;
    write_refinement_categories(data, block, *groups, &id)?;
    if groups.title_keywords {
        super::title_keywords::write_title_keywords_category(&state.info, &id, block);
    }
    Ok(id)
}
fn write_secondary_tail_categories(
    data: &BioStructureData,
    block: &mut CifBlock,
    groups: BioMmcifWriteParams,
    entry: &str,
    nontrivial_origx: bool,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   // _pdbx_struct_assembly* and _struct_biol are REMARK 300/350 in PDB
    // Gemmi❗❌:   if (groups.struct_biol && !st.meta.remark_300_detail.empty()) {
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_struct_biol.");
    // Gemmi❗❌:     span.set_pair("_struct_biol.id", "1");
    // Gemmi❗❌:     span.set_pair("_struct_biol.details", cif::quote(st.meta.remark_300_detail));
    // Gemmi❗❌:   }
    // Gemmi❗❌:   if (groups.assembly && !st.assemblies.empty())
    // Gemmi❗❌:     write_assemblies(st, block);
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.conn)
    // Gemmi❗❌:     write_struct_conn(st, block);
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.cis)  // _struct_mon_prot_cis
    // Gemmi❗❌:     write_cispeps(st, block);
    // Gemmi❗❌:
    // Gemmi❗❌:   // _pdbx_struct_mod_residue (MODRES)
    // Gemmi❗❌:   if (groups.modres && !st.mod_residues.empty()) {
    // Gemmi❗❌:     bool use_ccp4_mod_id = false;
    // Gemmi❗❌:     for (const ModRes& modres : st.mod_residues)
    // Gemmi❗❌:       if (!modres.mod_id.empty())
    // Gemmi❗❌:         use_ccp4_mod_id = true;
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_pdbx_struct_mod_residue.",
    // Gemmi❗❌:         {"id", "auth_asym_id", "auth_seq_id", "PDB_ins_code", "auth_comp_id",
    // Gemmi❗❌:          "label_comp_id", "parent_comp_id", "details"});
    // Gemmi❗❌:     if (use_ccp4_mod_id)
    // Gemmi❗❌:       loop.tags.push_back("_pdbx_struct_mod_residue.ccp4_mod_id");
    // Gemmi❗❌:     int counter = 0;
    // Gemmi❗❌:     for (const ModRes& modres : st.mod_residues) {
    // Gemmi❗❌:       loop.add_values({std::to_string(++counter),
    // Gemmi❗❌:                        qchain(modres.chain_name),
    // Gemmi❗❌:                        modres.res_id.seqid.num.str(),
    // Gemmi❗❌:                        pdbx_icode(modres.res_id),
    // Gemmi❗❌:                        string_or_dot(modres.res_id.name),
    // Gemmi❗❌:                        string_or_qmark(modres.res_id.name),
    // Gemmi❗❌:                        string_or_qmark(modres.parent_comp_id),
    // Gemmi❗❌:                        string_or_qmark(modres.details)});
    // Gemmi❗❌:       if (use_ccp4_mod_id)
    // Gemmi❗❌:         loop.values.push_back(string_or_qmark(modres.mod_id));
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   // _atom_sites (SCALE)
    // Gemmi❗❌:   if (groups.scale && (nontrivial_origx || st.cell.explicit_matrices)) {
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_atom_sites.");
    // Gemmi❗❌:     span.set_pair("_atom_sites.entry_id", id);
    // Gemmi❗❌:     std::string prefix = "_atom_sites.fract_transf_";
    // Gemmi❗❌:     for (int i = 0; i < 3; ++i) {
    // Gemmi❗❌:       std::string idx = "[" + std::to_string(i + 1) + "]";
    // Gemmi❗❌:       const auto& frac = st.cell.frac;
    // Gemmi❗❌:       std::string matrix_idx = prefix + "matrix";
    // Gemmi❗❌:       matrix_idx += idx;
    // Gemmi❗❌:       span.set_pair(matrix_idx + "[1]", to_str(frac.mat[i][0]));
    // Gemmi❗❌:       span.set_pair(matrix_idx + "[2]", to_str(frac.mat[i][1]));
    // Gemmi❗❌:       span.set_pair(matrix_idx + "[3]", to_str(frac.mat[i][2]));
    // Gemmi❗❌:       span.set_pair(cat(prefix, "vector", idx), to_str(frac.vec.at(i)));
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   // _atom_type
    // Gemmi❗❌:   if (groups.atom_type) {
    // Gemmi❗❌:     std::array<bool, (int)El::END> types{};
    // Gemmi❗❌:     for (const Model& model : st.models)
    // Gemmi❗❌:       for (const Chain& chain : model.chains)
    // Gemmi❗❌:         for (const Residue& res : chain.residues)
    // Gemmi❗❌:           for (const Atom& atom : res.atoms)
    // Gemmi❗❌:             types[atom.element.ordinal()] = true;
    // Gemmi❗❌:     cif::Loop& atom_type_loop = block.init_mmcif_loop("_atom_type.", {"symbol"});
    // Gemmi❗❌:     for (int i = 0; i < (int)El::END; ++i)
    // Gemmi❗❌:       if (types[i])
    // Gemmi❗❌:         atom_type_loop.add_row({Element((El)i).uname()});
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.entity_poly_seq) {
    // Gemmi❗❌:     cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.",
    // Gemmi❗❌:                                      {"entity_id", "num", "mon_id", "hetero"});
    // Gemmi❗❌:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       if (ent.entity_type == EntityType::Polymer) {
    // Gemmi❗❌:         // SEQRES from PDB doesn't record microheterogeneity.
    // Gemmi❗❌:         std::string hetero_no = ent.reflects_microhetero ? "n" : "?";
    // Gemmi❗❌:         for (size_t i = 0; i != ent.full_sequence.size(); ++i) {
    // Gemmi❗❌:           const std::string& mon_ids = ent.full_sequence[i];
    // Gemmi❗❌:           std::string num = std::to_string(i+1);
    // Gemmi❗❌:           size_t start = 0, end;
    // Gemmi❗❌:           while ((end = mon_ids.find(',', start)) != std::string::npos) {
    // Gemmi❗❌:             poly_loop.add_row({qchain(ent.name), num,
    // Gemmi❗❌:                                mon_ids.substr(start, end-start), "y"});
    // Gemmi❗❌:             start = end + 1;
    // Gemmi❗❌:           }
    // Gemmi❗❌:           poly_loop.add_row({qchain(ent.name), num, mon_ids.substr(start),
    // Gemmi❗❌:                              start == 0 ? hetero_no : "y"});
    // Gemmi❗❌:         }
    // Gemmi❗❌:       }
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   if (groups.atoms)
    // Gemmi❗❌:     add_cif_atoms(st, block, groups.group_pdb, groups.auth_all);
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let metadata = data.metadata();
    if groups.struct_biol && !metadata.remark_300_detail.is_empty() {
        block.set_pair_in_category(Some("_struct_biol."), "_struct_biol.id", "1".to_owned());
        block.set_pair_in_category(
            Some("_struct_biol."),
            "_struct_biol.details",
            quote_cif_value(&metadata.remark_300_detail),
        );
    }
    if groups.assembly && !data.assemblies().is_empty() {
        write_assemblies(data, block)?;
    }
    if groups.conn {
        write_struct_conn(data, block)?;
    }
    if groups.cis {
        write_cispeps(data, block)?;
    }
    if groups.modres && !data.mod_residues().is_empty() {
        let use_mod_id = data.mod_residues().iter().any(|m| !m.mod_id.is_empty());
        let mut tags = vec![
            "id",
            "auth_asym_id",
            "auth_seq_id",
            "PDB_ins_code",
            "auth_comp_id",
            "label_comp_id",
            "parent_comp_id",
            "details",
        ];
        if use_mod_id {
            tags.push("ccp4_mod_id");
        }
        let rows = data
            .mod_residues()
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let mut row = vec![
                    (i + 1).to_string(),
                    quote_cif_value(m.chain_name.as_str()),
                    m.res_id
                        .sequence_number()
                        .map_or_else(|| "?".to_owned(), |n| n.to_string()),
                    address_icode(&m.res_id),
                    string_or_dot(m.res_id.name().as_str()),
                    string_or_qmark(m.res_id.name().as_str()),
                    string_or_qmark(&m.parent_comp_id),
                    string_or_qmark(&m.details),
                ];
                if use_mod_id {
                    row.push(string_or_qmark(&m.mod_id));
                }
                row
            })
            .collect();
        add_mmcif_rows(block, "_pdbx_struct_mod_residue.", &tags, rows)?;
    }
    if groups.scale && (nontrivial_origx || data.crystal().is_some_and(|c| c.explicit_matrices())) {
        let identity = BioTransform::identity();
        let frac = data.crystal().map_or(&identity, |c| c.fractional());
        block.set_pair_in_category(
            Some("_atom_sites."),
            "_atom_sites.entry_id",
            entry.to_owned(),
        );
        for i in 0..3 {
            for j in 0..3 {
                block.set_pair_in_category(
                    Some("_atom_sites."),
                    &format!("_atom_sites.fract_transf_matrix[{}][{}]", i + 1, j + 1),
                    super::value::coordinate_text(frac.matrix()[i][j]),
                );
            }
            block.set_pair_in_category(
                Some("_atom_sites."),
                &format!("_atom_sites.fract_transf_vector[{}]", i + 1),
                super::value::coordinate_text(frac.translation()[i]),
            );
        }
    }
    if groups.atom_type {
        let mut elements = [false; 120];
        for atom in data.atoms() {
            elements[if atom.element() == cosmolkit_types::Element::H
                && atom.isotope_mass_number() == Some(2)
            {
                119
            } else {
                usize::from(atom.element().atomic_number())
            }] = true;
        }
        let rows = elements
            .into_iter()
            .enumerate()
            .filter(|(_, present)| *present)
            .map(|(i, _)| vec![cosmolkit_bio::GEMMI_ELEMENT_NAMES[i].to_ascii_uppercase()])
            .collect();
        add_mmcif_rows(block, "_atom_type.", &["symbol"], rows)?;
    }
    if groups.entity_poly_seq {
        let mut rows = Vec::new();
        for entity in data.entities() {
            if entity.kind() == EntityKind::Polymer {
                let hetero_no = if entity.reflects_microhetero() {
                    "n"
                } else {
                    "?"
                };
                for (index, monomers) in entity.full_sequence().iter().enumerate() {
                    let hetero = if monomers.contains(',') {
                        "y"
                    } else {
                        hetero_no
                    };
                    for monomer in monomers.split(',') {
                        rows.push(vec![
                            quote_cif_value(entity.source().source_entity_id()),
                            (index + 1).to_string(),
                            monomer.to_owned(),
                            hetero.to_owned(),
                        ]);
                    }
                }
            }
        }
        add_mmcif_rows(
            block,
            "_entity_poly_seq.",
            &["entity_id", "num", "mon_id", "hetero"],
            rows,
        )?;
    }
    if groups.atoms {
        super::atoms::add_cif_atoms(block, data, &groups)?;
    }
    write_tls_categories(data, block, groups)?;
    write_software_category(data, block, groups)
}
fn address_icode(address: &ResidueAddress) -> String {
    // Gemmi❗✔️: inline std::string pdbx_icode(const ResidueId& rid) {
    // Gemmi❗✔️:   return pdbx_icode(rid.seqid);
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.
    address
        .insertion_code()
        .map_or_else(|| "?".to_owned(), |b| char::from(b).to_string())
}

fn write_struct_ref_categories(
    data: &BioStructureData,
    block: &mut CifBlock,
    mapping: &std::collections::BTreeMap<&str, String>,
    entry: &str,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌:   if (groups.struct_ref) { // _struct_ref, _struct_ref_seq
    // Gemmi❗❌:     block.items.reserve(block.items.size() + 2); // avoid re-allocation
    // Gemmi❗❌:     cif::Loop& ref_loop = block.init_mmcif_loop("_struct_ref.",
    // Gemmi❗❌:                                   {"id", "entity_id", "db_name", "db_code",
    // Gemmi❗❌:                                    "pdbx_db_accession", "pdbx_db_isoform"});
    // Gemmi❗❌:     cif::Loop& seq_loop = block.init_mmcif_loop("_struct_ref_seq.", {
    // Gemmi❗❌:         "align_id", "ref_id", "pdbx_strand_id", "pdbx_PDB_id_code",
    // Gemmi❗❌:         "seq_align_beg", "seq_align_end", "pdbx_db_accession",
    // Gemmi❗❌:         "db_align_beg", "db_align_end",
    // Gemmi❗❌:         "pdbx_auth_seq_align_beg", "pdbx_seq_align_beg_ins_code",
    // Gemmi❗❌:         "pdbx_auth_seq_align_end", "pdbx_seq_align_end_ins_code"});
    // Gemmi❗❌:     int counter = 0;
    // Gemmi❗❌:     int counter2 = 0;
    // Gemmi❗❌:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       for (const Entity::DbRef& dbref : ent.dbrefs) {
    // Gemmi❗❌:         ref_loop.add_row({std::to_string(++counter),
    // Gemmi❗❌:                           qchain(ent.name),
    // Gemmi❗❌:                           string_or_dot(dbref.db_name),
    // Gemmi❗❌:                           string_or_dot(dbref.id_code),
    // Gemmi❗❌:                           string_or_qmark(dbref.accession_code),
    // Gemmi❗❌:                           string_or_qmark(dbref.isoform)});
    // Gemmi❗❌:         for (const std::string& subchain : ent.subchains) {
    // Gemmi❗❌:           auto strand_id = subs_to_strands.find(subchain);
    // Gemmi❗❌:           if (strand_id == subs_to_strands.end())
    // Gemmi❗❌:             continue;
    // Gemmi❗❌:           // DbRef::label_seq_begin/end (_struct_ref_seq.seq_align_beg/end) is
    // Gemmi❗❌:           // not filled in when reading PDB file, so we check it here.
    // Gemmi❗❌:           Residue::OptionalNum label_begin = dbref.label_seq_begin;
    // Gemmi❗❌:           Residue::OptionalNum label_end = dbref.label_seq_end;
    // Gemmi❗❌:           if (!label_begin || !label_end) {
    // Gemmi❗❌:             ConstResidueSpan span = st.models[0].get_subchain(subchain);
    // Gemmi❗❌:             try {
    // Gemmi❗❌:               label_begin = span.auth_seq_id_to_label(dbref.seq_begin);
    // Gemmi❗❌:               label_end = span.auth_seq_id_to_label(dbref.seq_end);
    // Gemmi❗❌:             } catch (const std::out_of_range&) {}
    // Gemmi❗❌:           }
    // Gemmi❗❌:           SeqId begin = dbref.seq_begin;
    // Gemmi❗❌:           SeqId end = dbref.seq_end;
    // Gemmi❗❌:           if (!begin.num || !end.num) {
    // Gemmi❗❌:             if (const Chain* chain = st.models[0].find_chain(strand_id->second))
    // Gemmi❗❌:               if (ConstResidueGroup polymer = chain->get_polymer()) {
    // Gemmi❗❌:                 begin = polymer.label_seq_id_to_auth(dbref.label_seq_begin);
    // Gemmi❗❌:                 end = polymer.label_seq_id_to_auth(dbref.label_seq_end);
    // Gemmi❗❌:               }
    // Gemmi❗❌:           }
    // Gemmi❗❌:           seq_loop.add_row({std::to_string(++counter2),
    // Gemmi❗❌:                             std::to_string(counter),
    // Gemmi❗❌:                             strand_id->second,  // pdbx_strand_id
    // Gemmi❗❌:                             id,
    // Gemmi❗❌:                             label_begin.str(),
    // Gemmi❗❌:                             label_end.str(),
    // Gemmi❗❌:                             string_or_qmark(dbref.accession_code),
    // Gemmi❗❌:                             dbref.db_begin.num.str(),
    // Gemmi❗❌:                             dbref.db_end.num.str(),
    // Gemmi❗❌:                             begin.num.str(),
    // Gemmi❗❌:                             pdbx_icode(begin),
    // Gemmi❗❌:                             end.num.str(),
    // Gemmi❗❌:                             pdbx_icode(end)});
    // Gemmi❗❌:         }
    // Gemmi❗❌:       }
    // Gemmi❗❌:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let mut refs = Vec::new();
    let mut seqs = Vec::new();
    for entity in data.entities() {
        for dbref in entity.dbrefs() {
            let ref_id = refs.len() + 1;
            refs.push(vec![
                ref_id.to_string(),
                quote_cif_value(entity.source().source_entity_id()),
                string_or_dot(&dbref.db_name),
                string_or_dot(&dbref.id_code),
                string_or_qmark(&dbref.accession_code),
                string_or_qmark(&dbref.isoform),
            ]);
            for sub in entity.subchains() {
                let Some(strand) = mapping.get(sub.as_str()) else {
                    continue;
                };
                let mut label_begin = dbref.label_seq_begin;
                let mut label_end = dbref.label_seq_end;
                if label_begin.is_none() || label_end.is_none() {
                    let residues = first_subchain(data, sub)?;
                    if !residues.is_empty() {
                        label_begin = auth_seq_id_to_label(&residues, Some(dbref.seq_begin))?;
                        label_end = auth_seq_id_to_label(&residues, Some(dbref.seq_end))?;
                    }
                }
                let mut begin = Some(dbref.seq_begin);
                let mut end = Some(dbref.seq_end);
                if dbref.seq_begin.seq_num() == i32::MIN || dbref.seq_end.seq_num() == i32::MIN {
                    if let Some(model) = data.models().first() {
                        if let Some(chain) =
                            model
                                .chain_span()
                                .slice(data.chains())?
                                .iter()
                                .find(|chain| {
                                    chain
                                        .source()
                                        .auth_chain_id()
                                        .is_some_and(|id| id.as_str() == strand)
                                })
                        {
                            let rows = chain.residue_span().slice(data.residues())?;
                            if let Some(start) = rows
                                .iter()
                                .position(|r| r.entity_kind() == EntityKind::Polymer)
                            {
                                let subchain = rows[start].source().subchain_id();
                                let count = rows[start..]
                                    .iter()
                                    .take_while(|r| {
                                        r.entity_kind() == EntityKind::Polymer
                                            && r.source().subchain_id() == subchain
                                    })
                                    .count();
                                let polymer = rows[start..start + count].iter().collect::<Vec<_>>();
                                begin = label_seq_id_to_auth(&polymer, dbref.label_seq_begin)?;
                                end = label_seq_id_to_auth(&polymer, dbref.label_seq_end)?;
                            }
                        }
                    }
                }
                seqs.push(vec![
                    (seqs.len() + 1).to_string(),
                    ref_id.to_string(),
                    strand.clone(),
                    entry.to_owned(),
                    label_begin.map_or_else(|| "?".to_owned(), |n| n.to_string()),
                    label_end.map_or_else(|| "?".to_owned(), |n| n.to_string()),
                    string_or_qmark(&dbref.accession_code),
                    seq_number_or_qmark(Some(dbref.db_begin)),
                    seq_number_or_qmark(Some(dbref.db_end)),
                    seq_number_or_qmark(begin),
                    pdbx_icode_value(begin),
                    seq_number_or_qmark(end),
                    pdbx_icode_value(end),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_struct_ref.",
        &[
            "id",
            "entity_id",
            "db_name",
            "db_code",
            "pdbx_db_accession",
            "pdbx_db_isoform",
        ],
        refs,
    )?;
    add_mmcif_rows(
        block,
        "_struct_ref_seq.",
        &[
            "align_id",
            "ref_id",
            "pdbx_strand_id",
            "pdbx_PDB_id_code",
            "seq_align_beg",
            "seq_align_end",
            "pdbx_db_accession",
            "db_align_beg",
            "db_align_end",
            "pdbx_auth_seq_align_beg",
            "pdbx_seq_align_beg_ins_code",
            "pdbx_auth_seq_align_end",
            "pdbx_seq_align_end_ins_code",
        ],
        seqs,
    )
}
fn first_subchain<'a>(
    data: &'a BioStructureData,
    name: &str,
) -> Result<Vec<&'a BioResidueRow>, BioMmcifWriteError> {
    // Gemmi❗❌:   ResidueSpan get_subchain(const std::string& sub_name) {
    // Gemmi❗❌:     for (Chain& chain : chains)
    // Gemmi❗❌:       if (ResidueSpan sub = chain.get_subchain(sub_name))
    // Gemmi❗❌:         return sub;
    // Gemmi❗❌:     return ResidueSpan();
    // Gemmi❗❌:   }
    // Gemmi❗❌:   ConstResidueSpan get_subchain(const std::string& sub_name) const {
    // Gemmi❗❌:     return const_cast<Model*>(this)->get_subchain(sub_name);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   ResidueSpan get_subchain(const std::string& s) {
    // Gemmi❗❌:     return get_residue_span([&](const Residue& r) { return r.subchain == s; });
    // Gemmi❗❌:   }
    // Gemmi❗❌:   ConstResidueSpan get_subchain(const std::string& s) const {
    // Gemmi❗❌:     return const_cast<Chain*>(this)->get_subchain(s);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   template<typename F> ResidueSpan get_residue_span(F&& func) {
    // Gemmi❗❌:     return whole().subspan(func);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   template<typename F> ConstResidueSpan get_residue_span(F&& func) const {
    // Gemmi❗❌:     return whole().subspan(func);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   ResidueSpan whole() {
    // Gemmi❗❌:     Residue* begin = residues.empty() ? nullptr : &residues[0];
    // Gemmi❗❌:     return ResidueSpan(residues, begin, residues.size());
    // Gemmi❗❌:   }
    // Gemmi❗❌:   ConstResidueSpan whole() const {
    // Gemmi❗❌:     const Residue* begin = residues.empty() ? nullptr : &residues[0];
    // Gemmi❗❌:     return ConstResidueSpan(begin, residues.size());
    // Gemmi❗❌:   }
    // Gemmi❗❌:   template<typename F, typename V=Item> Span<V> subspan(F&& func) {
    // Gemmi❗❌:     iterator group_begin = std::find_if(this->begin(), this->end(), func);
    // Gemmi❗❌:     iterator group_end = std::find_if_not(group_begin, this->end(), func);
    // Gemmi❗❌:     return Span<V>(&*group_begin, group_end - group_begin);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   template<typename F> Span<const value_type> subspan(F&& func) const {
    // Gemmi❗❌:     using V = const value_type;
    // Gemmi❗❌:     return const_cast<Span*>(this)->subspan<F, V>(std::forward<F>(func));
    // Gemmi❗❌:   }
    // Behavior: scan the first model and its chains in source order, then
    // return the first contiguous run with the requested subchain name. No
    // later chain or separated matching run is merged into this span.
    // Complexity: scans are O(N) as in the source helper chain. Gemmi returns
    // an O(1) borrowed span; Vec<&BioResidueRow> adds O(K) pointer copies and
    // allocation for K matching rows, with O(N) extra storage in the worst case.
    if let Some(model) = data.models().first() {
        for chain in model.chain_span().slice(data.chains())? {
            let rows = chain.residue_span().slice(data.residues())?;
            if let Some(start) = rows
                .iter()
                .position(|r| r.source().subchain_id().unwrap_or("") == name)
            {
                return Ok(rows[start..]
                    .iter()
                    .take_while(|r| r.source().subchain_id().unwrap_or("") == name)
                    .collect());
            }
        }
    }
    Ok(Vec::new())
}
fn write_assemblies(
    data: &BioStructureData,
    block: &mut CifBlock,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌: void write_assemblies(const Structure& st, cif::Block& block) {
    // Gemmi❗❌:   block.items.reserve(block.items.size() + 4); // avoid re-allocation
    // Gemmi❗❌:   cif::Loop& a_loop = block.init_mmcif_loop("_pdbx_struct_assembly.",
    // Gemmi❗❌:       {"id", "details", "method_details",
    // Gemmi❗❌:        "oligomeric_details", "oligomeric_count"});
    // Gemmi❗❌:   cif::Loop& prop_loop = block.init_mmcif_loop("_pdbx_struct_assembly_prop.",
    // Gemmi❗❌:       {"biol_id", "type", "value"});
    // Gemmi❗❌:   cif::Loop& gen_loop = block.init_mmcif_loop("_pdbx_struct_assembly_gen.",
    // Gemmi❗❌:       {"assembly_id", "oper_expression", "asym_id_list"});
    // Gemmi❗❌:   cif::Loop& oper_loop = block.init_mmcif_loop("_pdbx_struct_oper_list.",
    // Gemmi❗❌:       {"id", "type",
    // Gemmi❗❌:        "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
    // Gemmi❗❌:        "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
    // Gemmi❗❌:        "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
    // Gemmi❗❌:   std::vector<const Assembly::Operator*> distinct_oper;
    // Gemmi❗❌:   for (const Assembly& as : st.assemblies) {
    // Gemmi❗❌:     std::string how_defined = "?";
    // Gemmi❗❌:     if (as.author_determined && as.software_determined)
    // Gemmi❗❌:       how_defined = "author_and_software_defined_assembly";
    // Gemmi❗❌:     else if (as.author_determined)
    // Gemmi❗❌:       how_defined = "author_defined_assembly";
    // Gemmi❗❌:     else if (as.software_determined)
    // Gemmi❗❌:       how_defined = "software_defined_assembly";
    // Gemmi❗❌:     else if (as.special_kind == Assembly::SpecialKind::CompleteIcosahedral)
    // Gemmi❗❌:       how_defined = "'complete icosahedral assembly'";
    // Gemmi❗❌:     else if (as.special_kind == Assembly::SpecialKind::RepresentativeHelical)
    // Gemmi❗❌:       how_defined = "'representative helical assembly'";
    // Gemmi❗❌:     else if (as.special_kind == Assembly::SpecialKind::CompletePoint)
    // Gemmi❗❌:       how_defined = "'complete point assembly'";
    // Gemmi❗❌:     std::string oligomer = to_lower(as.oligomeric_details);
    // Gemmi❗❌:     int nmer = as.oligomeric_count != 0 ? as.oligomeric_count
    // Gemmi❗❌:                                         : xmeric_to_number(oligomer);
    // Gemmi❗❌:     // _pdbx_struct_assembly
    // Gemmi❗❌:     a_loop.add_row({as.name,
    // Gemmi❗❌:                     how_defined,
    // Gemmi❗❌:                     string_or_qmark(as.software_name),
    // Gemmi❗❌:                     string_or_qmark(oligomer),
    // Gemmi❗❌:                     nmer == 0 ? "?" : std::to_string(nmer)});
    // Gemmi❗❌:
    // Gemmi❗❌:     // _pdbx_struct_assembly_prop
    // Gemmi❗❌:     if (!std::isnan(as.absa))
    // Gemmi❗❌:       prop_loop.add_row({as.name, "'ABSA (A^2)'", to_str(as.absa)});
    // Gemmi❗❌:     if (!std::isnan(as.ssa))
    // Gemmi❗❌:       prop_loop.add_row({as.name, "'SSA (A^2)'", to_str(as.ssa)});
    // Gemmi❗❌:     if (!std::isnan(as.more))
    // Gemmi❗❌:       prop_loop.add_row({as.name, "MORE", to_str(as.more)});
    // Gemmi❗❌:
    // Gemmi❗❌:     // _pdbx_struct_assembly_gen and _pdbx_struct_oper_list
    // Gemmi❗❌:     for (const Assembly::Gen& gen : as.generators) {
    // Gemmi❗❌:       std::string subchain_str;
    // Gemmi❗❌:       for (const std::string& name : gen.subchains)
    // Gemmi❗❌:         string_append_sep(subchain_str, ',', name);
    // Gemmi❗❌:       if (subchain_str.empty()) // chain names to subchain names
    // Gemmi❗❌:         for (const Chain& chain : st.models[0].chains)
    // Gemmi❗❌:           if (in_vector(chain.name, gen.chains))
    // Gemmi❗❌:             for (const auto& sub : chain.subchains())
    // Gemmi❗❌:               string_append_sep(subchain_str, ',', sub.front().subchain);
    // Gemmi❗❌:       std::string oper_str;
    // Gemmi❗❌:       for (const Assembly::Operator& oper : gen.operators) {
    // Gemmi❗❌:         size_t k = 0;
    // Gemmi❗❌:         for (; k != distinct_oper.size(); ++k)
    // Gemmi❗❌:           if (distinct_oper[k]->transform.approx(oper.transform, 1e-9))
    // Gemmi❗❌:             break;
    // Gemmi❗❌:         string_append_sep(oper_str, ',', std::to_string(k+1));
    // Gemmi❗❌:         if (k != distinct_oper.size())
    // Gemmi❗❌:           continue;
    // Gemmi❗❌:         distinct_oper.emplace_back(&oper);
    // Gemmi❗❌:         oper_loop.values.emplace_back(std::to_string(k+1));
    // Gemmi❗❌:         if (!oper.type.empty()) {
    // Gemmi❗❌:           oper_loop.values.emplace_back(cif::quote(oper.type));
    // Gemmi❗❌:         } else if (oper.transform.is_identity()) {
    // Gemmi❗❌:           oper_loop.values.emplace_back("'identity operation'");
    // Gemmi❗❌:         } else if (as.author_determined || as.software_determined) {
    // Gemmi❗❌:           oper_loop.values.emplace_back("'crystal symmetry operation'");
    // Gemmi❗❌:         } else {
    // Gemmi❗❌:           oper_loop.values.emplace_back(".");
    // Gemmi❗❌:         }
    // Gemmi❗❌:         for (int i = 0; i < 3; ++i) {
    // Gemmi❗❌:           for (int j = 0; j < 3; ++j)
    // Gemmi❗❌:             oper_loop.values.emplace_back(to_str(oper.transform.mat[i][j]));
    // Gemmi❗❌:           oper_loop.values.emplace_back(to_str(oper.transform.vec.at(i)));
    // Gemmi❗❌:         }
    // Gemmi❗❌:       }
    // Gemmi❗❌:       gen_loop.add_row({as.name,
    // Gemmi❗❌:                         oper_str.empty() ? "." : oper_str,
    // Gemmi❗❌:                         subchain_str.empty() ? "?" : subchain_str});
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let mut a_rows = Vec::new();
    let mut prop_rows = Vec::new();
    let mut gen_rows = Vec::new();
    let mut oper_rows = Vec::new();
    let mut distinct = Vec::<&BioAssemblyOperator>::new();
    for assembly in data.assemblies() {
        let how = if assembly.author_determined && assembly.software_determined {
            "author_and_software_defined_assembly"
        } else if assembly.author_determined {
            "author_defined_assembly"
        } else if assembly.software_determined {
            "software_defined_assembly"
        } else {
            match assembly.special_kind {
                BioAssemblySpecialKind::CompleteIcosahedral => "'complete icosahedral assembly'",
                BioAssemblySpecialKind::RepresentativeHelical => {
                    "'representative helical assembly'"
                }
                BioAssemblySpecialKind::CompletePoint => "'complete point assembly'",
                BioAssemblySpecialKind::NotApplicable => "?",
            }
        };
        let oligomer = assembly.oligomeric_details.to_ascii_lowercase();
        let nmer = if assembly.oligomeric_count != 0 {
            assembly.oligomeric_count
        } else {
            xmeric_to_number(&oligomer)
        };
        a_rows.push(vec![
            assembly.name.clone(),
            how.to_owned(),
            string_or_qmark(&assembly.software_name),
            string_or_qmark(&oligomer),
            if nmer == 0 {
                "?".to_owned()
            } else {
                nmer.to_string()
            },
        ]);
        for (kind, value) in [
            ("'ABSA (A^2)'", assembly.buried_surface_area),
            ("'SSA (A^2)'", assembly.surface_area),
            ("MORE", assembly.solvent_free_energy_change),
        ] {
            if !value.is_nan() {
                prop_rows.push(vec![
                    assembly.name.clone(),
                    kind.to_owned(),
                    super::value::coordinate_text(value),
                ]);
            }
        }
        for generator in &assembly.generators {
            let mut subs = generator.subchains.join(",");
            if subs.is_empty() {
                if let Some(model) = data.models().first() {
                    for chain in model.chain_span().slice(data.chains())? {
                        if chain
                            .source()
                            .auth_chain_id()
                            .is_some_and(|id| generator.chains.iter().any(|s| s == id.as_str()))
                        {
                            let mut previous = None;
                            for residue in chain.residue_span().slice(data.residues())? {
                                let sub = residue.source().subchain_id().unwrap_or("");
                                if previous == Some(sub) {
                                    continue;
                                }
                                previous = Some(sub);
                                if !subs.is_empty() {
                                    subs.push(',');
                                }
                                subs.push_str(sub);
                            }
                        }
                    }
                }
            }
            let mut op_ids = Vec::new();
            for operator in &generator.operators {
                let k = distinct
                    .iter()
                    .position(|candidate| {
                        transform_approx(&candidate.transform, &operator.transform)
                    })
                    .unwrap_or(distinct.len());
                op_ids.push((k + 1).to_string());
                if k != distinct.len() {
                    continue;
                }
                distinct.push(operator);
                let kind = if let Some(kind) =
                    operator.operator_type.as_deref().filter(|s| !s.is_empty())
                {
                    quote_cif_value(kind)
                } else if transform_is_exact_identity(&operator.transform) {
                    "'identity operation'".to_owned()
                } else if assembly.author_determined || assembly.software_determined {
                    "'crystal symmetry operation'".to_owned()
                } else {
                    ".".to_owned()
                };
                let mut row = vec![(k + 1).to_string(), kind];
                for i in 0..3 {
                    for j in 0..3 {
                        row.push(super::value::coordinate_text(
                            operator.transform.matrix()[i][j],
                        ));
                    }
                    row.push(super::value::coordinate_text(
                        operator.transform.translation()[i],
                    ));
                }
                oper_rows.push(row);
            }
            gen_rows.push(vec![
                assembly.name.clone(),
                if op_ids.is_empty() {
                    ".".to_owned()
                } else {
                    op_ids.join(",")
                },
                if subs.is_empty() {
                    "?".to_owned()
                } else {
                    subs
                },
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_struct_assembly.",
        &[
            "id",
            "details",
            "method_details",
            "oligomeric_details",
            "oligomeric_count",
        ],
        a_rows,
    )?;
    add_mmcif_rows(
        block,
        "_pdbx_struct_assembly_prop.",
        &["biol_id", "type", "value"],
        prop_rows,
    )?;
    add_mmcif_rows(
        block,
        "_pdbx_struct_assembly_gen.",
        &["assembly_id", "oper_expression", "asym_id_list"],
        gen_rows,
    )?;
    add_mmcif_rows(
        block,
        "_pdbx_struct_oper_list.",
        &[
            "id",
            "type",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "vector[1]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "vector[2]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[3]",
        ],
        oper_rows,
    )
}
fn xmeric_to_number(oligomer: &str) -> i32 {
    // Gemmi❗✔️: int xmeric_to_number(const std::string& oligomeric) {
    // Gemmi❗✔️:   static const char names[20][10] = {
    // Gemmi❗✔️:     "mono", "di", "tri", "tetra", "penta",
    // Gemmi❗✔️:     "hexa", "hepta", "octa", "nona", "deca",
    // Gemmi❗✔️:     "undeca", "dodeca", "trideca", "tetradeca", "pentadeca",
    // Gemmi❗✔️:     "hexadeca", "heptadeca", "octadeca", "nonadeca", "eicosa"
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   size_t len = oligomeric.length();
    // Gemmi❗✔️:   const char* p = oligomeric.c_str();
    // Gemmi❗✔️:   for (int i = 0; i != 20; ++i)
    // Gemmi❗✔️:     if (len == std::strlen(names[i]) + 5 && strncmp(p, names[i], len-5) == 0)
    // Gemmi❗✔️:       return i + 1;
    // Gemmi❗✔️:   return no_sign_atoi(p);
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    // Gemmi❗✔️: inline int no_sign_atoi(const char* p, const char** endptr=nullptr) {
    // Gemmi❗✔️:   int n = 0;
    // Gemmi❗✔️:   while (is_space(*p))
    // Gemmi❗✔️:     ++p;
    // Gemmi❗✔️:   for (; is_digit(*p); ++p)
    // Gemmi❗✔️:     n = n * 10 + (*p - '0');
    // Gemmi❗✔️:   if (endptr)
    // Gemmi❗✔️:     *endptr = p;
    // Gemmi❗✔️:   return n;
    // Gemmi❗✔️: }
    const NAMES: [&str; 20] = [
        "mono",
        "di",
        "tri",
        "tetra",
        "penta",
        "hexa",
        "hepta",
        "octa",
        "nona",
        "deca",
        "undeca",
        "dodeca",
        "trideca",
        "tetradeca",
        "pentadeca",
        "hexadeca",
        "heptadeca",
        "octadeca",
        "nonadeca",
        "eicosa",
    ];
    for (i, prefix) in NAMES.iter().enumerate() {
        if oligomer.len() == prefix.len() + 5 && oligomer.starts_with(prefix) {
            return i as i32 + 1;
        }
    }
    let mut n = 0_i32;
    for digit in oligomer
        .trim_start_matches(|c: char| matches!(c, '\u{9}'..='\u{d}' | ' '))
        .bytes()
        .take_while(u8::is_ascii_digit)
    {
        n = n.wrapping_mul(10).wrapping_add(i32::from(digit - b'0'));
    }
    n
}
fn connection_type_text(kind: BioConnectionKind) -> &'static str {
    // Gemmi❗✔️: inline const char* connection_type_to_string(Connection::Type t) {
    // Gemmi❗✔️:   static constexpr const char* type_ids[] = {
    // Gemmi❗✔️:     "covale", "disulf", "hydrog", "metalc", "."
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   return type_ids[t];
    // Gemmi❗✔️: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.
    match kind {
        BioConnectionKind::Covale => "covale",
        BioConnectionKind::Disulf => "disulf",
        BioConnectionKind::Hydrog => "hydrog",
        BioConnectionKind::MetalC => "metalc",
        BioConnectionKind::Unknown => ".",
    }
}
fn connection_values(
    data: &BioStructureData,
    residue: BioResidueId,
    atom: Option<BioAtomId>,
    address: &AtomAddress,
) -> Vec<String> {
    // Gemmi❗❌:     v.emplace_back(subchain_or_dot(*cra1.residue));       // ptnr1_label_asym_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->name);                   // ptnr1_label_comp_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->label_seq.str('.'));     // ptnr1_label_seq_id
    // Gemmi❗❌:     v.emplace_back(at1 ? cif::quote(at1->name) : "?");    // ptnr1_label_atom_id
    // Gemmi❗❌:     v.emplace_back(1, at1 ? at1->altloc_or('?') : '?');   // pdbx_ptnr1_label_alt_id
    // Gemmi❗❌:     v.emplace_back(qchain(con.partner1.chain_name));      // ptnr1_auth_asym_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->seqid.num.str());        // ptnr1_auth_seq_id
    // Gemmi❗❌:     v.emplace_back(pdbx_icode(con.partner1.res_id));      // ptnr1_PDB_ins_code
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let residue = &data.residues()[residue.index()];
    let atom = atom.map(|id| &data.atoms()[id.index()]);
    vec![
        string_or_dot(residue.source().subchain_id().unwrap_or("")),
        residue_name_logical_view(&residue.name(), data.input_format()).to_owned(),
        residue
            .source()
            .label_seq_id()
            .map_or_else(|| ".".to_owned(), |n| n.to_string()),
        atom.map_or_else(
            || "?".to_owned(),
            |a| quote_cif_value(atom_name_logical_view(&a.name(), data.input_format())),
        ),
        atom.and_then(|a| a.altloc())
            .map_or_else(|| "?".to_owned(), |alt| char::from(alt.value()).to_string()),
        quote_cif_value(address.chain_name().as_str()),
        seq_number_or_qmark(residue.source().seq_id()),
        address_icode(&address.residue()),
    ]
}
fn symmetry_code(image: &BioNearestImage) -> String {
    // Gemmi❗✔️:   std::string symmetry_code(bool underscore) const {
    // Gemmi❗✔️:     std::string s = std::to_string(sym_idx + 1);
    // Gemmi❗✔️:     if (underscore)
    // Gemmi❗✔️:       s += '_';
    // Gemmi❗✔️:     if (unsigned(5 + pbc_shift[0]) <= 9 &&
    // Gemmi❗✔️:         unsigned(5 + pbc_shift[1]) <= 9 &&
    // Gemmi❗✔️:         unsigned(5 + pbc_shift[2]) <= 9) {  // normal, quick path
    // Gemmi❗✔️:       for (int shift : pbc_shift)
    // Gemmi❗✔️:         s += char('5' + shift);
    // Gemmi❗✔️:     } else {                                // problematic, non-standard path
    // Gemmi❗✔️:       for (int i = 0; i < 3; ++i) {
    // Gemmi❗✔️:         if (i != 0 && underscore)
    // Gemmi❗✔️:           s += '_';
    // Gemmi❗✔️:         s += std::to_string(5 + pbc_shift[i]);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     return s;
    // Gemmi❗✔️:   }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: scalar/source-shaped traversal, without extra row buffering.

    let mut s = (image.sym_idx() + 1).to_string();
    s.push('_');
    if image
        .pbc_shift()
        .iter()
        .all(|shift| (0..=9).contains(&(5 + shift)))
    {
        for shift in image.pbc_shift() {
            s.push(char::from((b'5' as i32 + shift) as u8));
        }
    } else {
        for (i, shift) in image.pbc_shift().iter().enumerate() {
            if i != 0 {
                s.push('_');
            }
            s.push_str(&(5 + shift).to_string());
        }
    }
    s
}
fn default_crystal() -> BioCrystalInfo {
    // Gemmi❗✔️:   double a = 1.0, b = 1.0, c = 1.0;
    // Gemmi❗✔️:   double alpha = 90.0, beta = 90.0, gamma = 90.0;

    BioCrystalInfo::new(
        BioCrystalCell::default(),
        None,
        None,
        BioTransform::identity(),
        BioTransform::identity(),
        false,
        0,
        Vec::new(),
    )
}
fn write_struct_conn(
    data: &BioStructureData,
    block: &mut CifBlock,
) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌: void write_struct_conn(const Structure& st, cif::Block& block) {
    // Gemmi❗❌:   // example:
    // Gemmi❗❌:   // disulf1 disulf A CYS 3  SG ? 3 ? 1_555 A CYS 18 SG ? 18 ?  1_555 ? 2.045
    // Gemmi❗❌:   std::array<bool,(int)Connection::Type::Unknown+1> type_ids{};
    // Gemmi❗❌:   bool use_ccp4_link_id = false;
    // Gemmi❗❌:   for (const Connection& con : st.connections)
    // Gemmi❗❌:     if (!con.link_id.empty())
    // Gemmi❗❌:       use_ccp4_link_id = true;
    // Gemmi❗❌:   cif::Loop& conn_loop = block.init_mmcif_loop("_struct_conn.",
    // Gemmi❗❌:       {"id", "conn_type_id",
    // Gemmi❗❌:        "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id",
    // Gemmi❗❌:        "ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id", "ptnr1_auth_asym_id",
    // Gemmi❗❌:        "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_symmetry",
    // Gemmi❗❌:        "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id",
    // Gemmi❗❌:        "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id", "ptnr2_auth_asym_id",
    // Gemmi❗❌:        "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code", "ptnr2_symmetry",
    // Gemmi❗❌:        "details", "pdbx_dist_value"});
    // Gemmi❗❌:   if (use_ccp4_link_id)
    // Gemmi❗❌:     conn_loop.tags.push_back("_struct_conn.ccp4_link_id");
    // Gemmi❗❌:   for (const Connection& con : st.connections) {
    // Gemmi❗❌:     const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
    // Gemmi❗❌:     const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
    // Gemmi❗❌:     if (!cra1.residue || !cra2.residue)
    // Gemmi❗❌:       continue;
    // Gemmi❗❌:     const Atom* at1 = cra1.atom;
    // Gemmi❗❌:     const Atom* at2 = cra2.atom;
    // Gemmi❗❌:     std::string im_pdb_symbol = "?", im_dist_str = "?";
    // Gemmi❗❌:     if (at1 && at2) {
    // Gemmi❗❌:       NearestImage im = st.cell.find_nearest_image(at1->pos, at2->pos, con.asu);
    // Gemmi❗❌:       im_pdb_symbol = im.symmetry_code(true);
    // Gemmi❗❌:       im_dist_str = to_str_prec<4>(im.dist());
    // Gemmi❗❌:     }
    // Gemmi❗❌:     auto& v = conn_loop.values;
    // Gemmi❗❌:     v.emplace_back(string_or_qmark(con.name));            // id
    // Gemmi❗❌:     v.emplace_back(connection_type_to_string(con.type));  // conn_type_id
    // Gemmi❗❌:     v.emplace_back(subchain_or_dot(*cra1.residue));       // ptnr1_label_asym_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->name);                   // ptnr1_label_comp_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->label_seq.str('.'));     // ptnr1_label_seq_id
    // Gemmi❗❌:     v.emplace_back(at1 ? cif::quote(at1->name) : "?");    // ptnr1_label_atom_id
    // Gemmi❗❌:     v.emplace_back(1, at1 ? at1->altloc_or('?') : '?');   // pdbx_ptnr1_label_alt_id
    // Gemmi❗❌:     v.emplace_back(qchain(con.partner1.chain_name));      // ptnr1_auth_asym_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->seqid.num.str());        // ptnr1_auth_seq_id
    // Gemmi❗❌:     v.emplace_back(pdbx_icode(con.partner1.res_id));      // ptnr1_PDB_ins_code
    // Gemmi❗❌:     v.emplace_back("1_555");                              // ptnr1_symmetry
    // Gemmi❗❌:     v.emplace_back(subchain_or_dot(*cra2.residue));       // ptnr2_label_asym_id
    // Gemmi❗❌:     v.emplace_back(cra2.residue->name);                   // ptnr2_label_comp_id
    // Gemmi❗❌:     v.emplace_back(cra2.residue->label_seq.str('.'));     // ptnr2_label_seq_id
    // Gemmi❗❌:     v.emplace_back(at2 ? cif::quote(at2->name) : "?");    // ptnr2_label_atom_id
    // Gemmi❗❌:     v.emplace_back(1, at2 ? at2->altloc_or('?') : '?');   // pdbx_ptnr2_label_alt_id
    // Gemmi❗❌:     v.emplace_back(qchain(con.partner2.chain_name));      // ptnr2_auth_asym_id
    // Gemmi❗❌:     v.emplace_back(cra2.residue->seqid.num.str());        // ptnr2_auth_seq_id
    // Gemmi❗❌:     v.emplace_back(pdbx_icode(con.partner2.res_id));      // ptnr2_PDB_ins_code
    // Gemmi❗❌:     v.emplace_back(im_pdb_symbol);                        // ptnr2_symmetry
    // Gemmi❗❌:     v.emplace_back("?");                                  // details
    // Gemmi❗❌:     v.emplace_back(im_dist_str);                          // pdbx_dist_value
    // Gemmi❗❌:     if (use_ccp4_link_id)
    // Gemmi❗❌:       v.emplace_back(string_or_qmark(con.link_id));       // ccp4_link_id
    // Gemmi❗❌:     type_ids[int(con.type)] = true;
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   cif::Loop& type_loop = block.init_mmcif_loop("_struct_conn_type.", {"id"});
    // Gemmi❗❌:   for (int i = 0; i < (int)type_ids.size() - 1; ++i)
    // Gemmi❗❌:     if (type_ids[i])
    // Gemmi❗❌:       type_loop.add_row({connection_type_to_string((Connection::Type)i)});
    // Gemmi❗❌: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let link_id = data.connections().iter().any(|c| !c.link_id.is_empty());
    let mut tags = vec![
        "id",
        "conn_type_id",
        "ptnr1_label_asym_id",
        "ptnr1_label_comp_id",
        "ptnr1_label_seq_id",
        "ptnr1_label_atom_id",
        "pdbx_ptnr1_label_alt_id",
        "ptnr1_auth_asym_id",
        "ptnr1_auth_seq_id",
        "pdbx_ptnr1_PDB_ins_code",
        "ptnr1_symmetry",
        "ptnr2_label_asym_id",
        "ptnr2_label_comp_id",
        "ptnr2_label_seq_id",
        "ptnr2_label_atom_id",
        "pdbx_ptnr2_label_alt_id",
        "ptnr2_auth_asym_id",
        "ptnr2_auth_seq_id",
        "pdbx_ptnr2_PDB_ins_code",
        "ptnr2_symmetry",
        "details",
        "pdbx_dist_value",
    ];
    if link_id {
        tags.push("ccp4_link_id");
    }
    let mut rows = Vec::new();
    let mut types = [false; 5];
    let default = default_crystal();
    let crystal = data.crystal().unwrap_or(&default);
    for connection in data.connections() {
        let Some((_, left_res, left_at)) =
            data.find_cra(BioModelId::new(0), &connection.partner1, true)?
        else {
            continue;
        };
        let Some((_, right_res, right_at)) =
            data.find_cra(BioModelId::new(0), &connection.partner2, true)?
        else {
            continue;
        };
        let (symmetry, distance) = if let (Some(left), Some(right)) = (left_at, right_at) {
            let image = cosmolkit_bio::find_nearest_image(
                crystal,
                data.coordinates().positions()[left.index()],
                data.coordinates().positions()[right.index()],
                connection.asu,
            );
            (
                symmetry_code(&image),
                crate::cif::format_cif_f64_precision::<4>(image.dist_sq().sqrt()),
            )
        } else {
            ("?".to_owned(), "?".to_owned())
        };
        let mut row = vec![
            string_or_qmark(&connection.name),
            connection_type_text(connection.kind).to_owned(),
        ];
        row.extend(connection_values(
            data,
            left_res,
            left_at,
            &connection.partner1,
        ));
        row.push("1_555".to_owned());
        row.extend(connection_values(
            data,
            right_res,
            right_at,
            &connection.partner2,
        ));
        row.extend([symmetry, "?".to_owned(), distance]);
        if link_id {
            row.push(string_or_qmark(&connection.link_id));
        }
        rows.push(row);
        types[connection.kind as usize] = true;
    }
    add_mmcif_rows(block, "_struct_conn.", &tags, rows)?;
    add_mmcif_rows(
        block,
        "_struct_conn_type.",
        &["id"],
        [
            BioConnectionKind::Covale,
            BioConnectionKind::Disulf,
            BioConnectionKind::Hydrog,
            BioConnectionKind::MetalC,
        ]
        .into_iter()
        .filter(|kind| types[*kind as usize])
        .map(|kind| vec![connection_type_text(kind).to_owned()])
        .collect(),
    )
}
fn write_cispeps(data: &BioStructureData, block: &mut CifBlock) -> Result<(), BioMmcifWriteError> {
    // Gemmi❗❌: void write_cispeps(const Structure& st, cif::Block& block) {
    // Gemmi❗❌:   cif::Loop* prot_cis_loop = nullptr;
    // Gemmi❗❌:   int pdbx_id = 0;
    // Gemmi❗❌:   for (const CisPep& cispep : st.cispeps) {
    // Gemmi❗❌:     const Model* model = &st.models[0];
    // Gemmi❗❌:     if (st.models.size() > 1) {
    // Gemmi❗❌:       model = st.find_model(cispep.model_num);
    // Gemmi❗❌:       if (!model)
    // Gemmi❗❌:         continue;
    // Gemmi❗❌:     }
    // Gemmi❗❌:     const_CRA cra1 = model->find_cra(cispep.partner_c, true);
    // Gemmi❗❌:     const_CRA cra2 = model->find_cra(cispep.partner_n, true);
    // Gemmi❗❌:     if (!cra1.residue || !cra2.residue)
    // Gemmi❗❌:       continue;
    // Gemmi❗❌:     if (!prot_cis_loop)
    // Gemmi❗❌:       prot_cis_loop = &block.init_mmcif_loop("_struct_mon_prot_cis.",
    // Gemmi❗❌:           {"pdbx_id", "pdbx_PDB_model_num",
    // Gemmi❗❌:            "label_asym_id", "label_seq_id", "label_comp_id",
    // Gemmi❗❌:            "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
    // Gemmi❗❌:            "pdbx_label_asym_id_2", "pdbx_label_seq_id_2", "pdbx_label_comp_id_2",
    // Gemmi❗❌:            "pdbx_auth_asym_id_2", "pdbx_auth_seq_id_2", "pdbx_PDB_ins_code_2",
    // Gemmi❗❌:            "label_alt_id", "pdbx_omega_angle"});
    // Gemmi❗❌:     auto& v = prot_cis_loop->values;
    // Gemmi❗❌:     v.emplace_back(std::to_string(++pdbx_id));            // pdbx_id
    // Gemmi❗❌:     v.emplace_back(std::to_string(model->num));           // pdbx_PDB_model_num
    // Gemmi❗❌:     v.emplace_back(subchain_or_dot(*cra1.residue));       // label_asym_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->label_seq.str('.'));     // label_seq_id
    // Gemmi❗❌:     v.emplace_back(cra1.residue->name);                   // label_comp_id
    // Gemmi❗❌:     v.emplace_back(qchain(cispep.partner_c.chain_name));  // auth_asym_id
    // Gemmi❗❌:     v.emplace_back(cispep.partner_c.res_id.seqid.num.str()); // auth_seq_id
    // Gemmi❗❌:     v.emplace_back(pdbx_icode(cispep.partner_c.res_id));  // pdbx_PDB_ins_code
    // Gemmi❗❌:     v.emplace_back(subchain_or_dot(*cra2.residue));       // pdbx_label_asym_id_2
    // Gemmi❗❌:     v.emplace_back(cra2.residue->label_seq.str('.'));     // pdbx_label_seq_id_2
    // Gemmi❗❌:     v.emplace_back(cra2.residue->name);                   // pdbx_label_comp_id_2
    // Gemmi❗❌:     v.emplace_back(qchain(cispep.partner_n.chain_name));  // pdbx_auth_asym_id_2
    // Gemmi❗❌:     v.emplace_back(cispep.partner_n.res_id.seqid.num.str()); // pdbx_auth_seq_id_2
    // Gemmi❗❌:     v.emplace_back(pdbx_icode(cispep.partner_n.res_id));  // pdbx_PDB_ins_code_2
    // Gemmi❗❌:     v.emplace_back(1, cispep.only_altloc ? cispep.only_altloc : '.');
    // Gemmi❗❌:     v.emplace_back(number_or_qmark(cispep.reported_angle));
    // Gemmi❗❌:   }
    // Gemmi❗❌: }
    // Behavior: source-shaped implementation; native full-profile verification is retained separately.
    // Cost: temporary row vectors or repeated category scans add allocations/work versus direct source append.

    let mut rows = Vec::new();
    for cis in data.cispeps() {
        let model = if data.models().len() > 1 {
            let Some(index) = data
                .models()
                .iter()
                .position(|m| m.source_model_number() == Some(cis.model_num))
            else {
                continue;
            };
            index
        } else {
            0
        };
        let Some((_, left_res, left_at)) =
            data.find_cra(BioModelId::new(model as u32), &cis.partner_c, true)?
        else {
            continue;
        };
        let Some((_, right_res, right_at)) =
            data.find_cra(BioModelId::new(model as u32), &cis.partner_n, true)?
        else {
            continue;
        };
        let left = connection_values(data, left_res, left_at, &cis.partner_c);
        let right = connection_values(data, right_res, right_at, &cis.partner_n);
        rows.push(vec![
            (rows.len() + 1).to_string(),
            data.models()[model]
                .source_model_number()
                .unwrap_or_default()
                .to_string(),
            left[0].clone(),
            left[2].clone(),
            left[1].clone(),
            left[5].clone(),
            cis.partner_c
                .residue()
                .sequence_number()
                .map_or_else(|| "?".to_owned(), |n| n.to_string()),
            left[7].clone(),
            right[0].clone(),
            right[2].clone(),
            right[1].clone(),
            right[5].clone(),
            cis.partner_n
                .residue()
                .sequence_number()
                .map_or_else(|| "?".to_owned(), |n| n.to_string()),
            right[7].clone(),
            char::from(if cis.only_altloc == 0 {
                b'.'
            } else {
                cis.only_altloc
            })
            .to_string(),
            number_or_qmark(Some(cis.reported_angle)),
        ]);
    }
    if rows.is_empty() {
        return Ok(());
    }
    add_mmcif_rows(
        block,
        "_struct_mon_prot_cis.",
        &[
            "pdbx_id",
            "pdbx_PDB_model_num",
            "label_asym_id",
            "label_seq_id",
            "label_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "pdbx_PDB_ins_code",
            "pdbx_label_asym_id_2",
            "pdbx_label_seq_id_2",
            "pdbx_label_comp_id_2",
            "pdbx_auth_asym_id_2",
            "pdbx_auth_seq_id_2",
            "pdbx_PDB_ins_code_2",
            "label_alt_id",
            "pdbx_omega_angle",
        ],
        rows,
    )
}
