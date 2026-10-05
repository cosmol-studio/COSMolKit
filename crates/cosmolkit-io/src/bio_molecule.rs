//! BIO to detached molecular blocks; reuses the existing PDB chemistry owner.
use cosmolkit_bio::{BioStructureData, BioStructureError};
use cosmolkit_model::{
    Atom, AtomId, AtomPdbResidueInfo, AtomSpec, Bond, Conformer3D, CoordinateBlock,
    MoleculeProperties, TopologyBlock,
};
use std::collections::HashMap;

/// Canonical detached structural conversion controls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BioMoleculeParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub flavor: u32,
    pub proximity_bonding: bool,
}
impl Default for BioMoleculeParams {
    fn default() -> Self {
        // RDKit❗✔️: struct RDKIT_FILEPARSERS_EXPORT PDBParserParams {
        // RDKit❗✔️:   bool sanitize = true; /**< sanitize the molecule after building it */
        // RDKit❗✔️:   bool removeHs = true; /**< remove Hs after constructing the molecule */
        // RDKit❗✔️:   bool proximityBonding = true; /**< if set to true, proximity bonding will be
        // RDKit❗✔️:                                    performed */
        // RDKit❗✔️:   unsigned int flavor = 0;      /**< flavor to use */
        // RDKit❗✔️: };
        Self {
            sanitize: true,
            remove_hs: true,
            flavor: 0,
            proximity_bonding: true,
        }
    }
}
#[derive(Debug, thiserror::Error)]
pub enum BioMoleculeConversionError {
    #[error(transparent)]
    Structure(#[from] BioStructureError),
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Chemistry(#[from] crate::PdbPostprocessError),
    #[error(transparent)]
    Sanitize(#[from] cosmolkit_core::SanitizeError),
    #[error(transparent)]
    Hydrogens(#[from] cosmolkit_core::HydrogenError),
    #[error(transparent)]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    Stereo(#[from] cosmolkit_core::StereoError),
}

/// Convert all hierarchy rows in encounter order to one molecular conformer.
/// Input parsing has already established the element and isotope values.
/// Missing original PDB text on manually built or CIF rows supplies no lexical
/// pseudo-atom predicate; actual PDB rows retain the original 24 bytes.
pub fn bio_structure_to_molecule_parts(
    data: &BioStructureData,
    params: &BioMoleculeParams,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), BioMoleculeConversionError> {
    // RDKit❗✔️:
    // RDKit❗✔️:   if ((flavor & 1) == 0) {
    // RDKit❗✔️:     // Ignore alternate locations of atoms.
    // RDKit❗✔️:     if (len >= 17 && ptr[16] != ' ' && ptr[16] != 'A' && ptr[16] != '1') {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // Ignore XPLOR pseudo atoms
    // RDKit❗✔️:     if (len >= 54 && !memcmp(ptr + 30, "9999.0009999.0009999.000", 24)) {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // Ignore NMR pseudo atoms
    // RDKit❗✔️:     if (ptr[12] == ' ' && ptr[13] == 'Q') {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // Ignore PDB dummy residues
    // RDKit❗✔️:     if (len >= 20 && !memcmp(ptr + 18, "DUM", 3)) {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    data.validate()?;
    let mut atoms = Vec::with_capacity(data.atoms().len());
    let mut positions = Vec::with_capacity(data.atoms().len());
    let mut serials = HashMap::new();
    for (index, row) in data.atoms().iter().enumerate() {
        let residue = &data.residues()[row.residue_id().index()];
        let chain = &data.chains()[residue.chain_id().index()];
        let atom_name = row.name();
        let name = atom_name.as_bytes();
        if params.flavor & 1 == 0
            && (row
                .altloc()
                .is_some_and(|a| a.value() != b'A' && a.value() != b'1')
                || row
                    .pdb_coordinate_text()
                    .is_some_and(|text| text == b"9999.0009999.0009999.000")
                || (name[0] == b' ' && name[1] == b'Q')
                || residue.name().as_str() == "DUM")
        {
            continue;
        }
        let id = AtomId::new(atoms.len());
        let serial = row.source().serial().map_or(0, |x| x.value());
        let seq = residue.source().seq_id();
        // RDKit❗✔️:   AtomPDBResidueInfo *info = new AtomPDBResidueInfo(tmp, serialno);
        // RDKit❗✔️:   atom->setMonomerInfo(info);
        // RDKit❗✔️:
        // RDKit❗✔️:   if (len >= 20) {
        // RDKit❗✔️:     tmp = std::string(ptr + 17, 3);
        // RDKit❗✔️:     // boost::trim(tmp);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     tmp = "UNL";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setResidueName(tmp);
        // RDKit❗✔️:   if (ptr[0] == 'H') {
        // RDKit❗✔️:     info->setIsHeteroAtom(true);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (len >= 17) {
        // RDKit❗✔️:     tmp = std::string(ptr + 16, 1);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     tmp = " ";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setAltLoc(tmp);
        // RDKit❗✔️:   if (len >= 22) {
        // RDKit❗✔️:     tmp = std::string(ptr + 21, 1);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     tmp = " ";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setChainId(tmp);
        // RDKit❗✔️:   if (len >= 27) {
        // RDKit❗✔️:     tmp = std::string(ptr + 26, 1);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     tmp = " ";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setInsertionCode(tmp);
        // RDKit❗✔️:
        // RDKit❗✔️:   int resno = 1;
        // RDKit❗✔️:   if (len >= 26) {
        // RDKit❗✔️:     try {
        // RDKit❗✔️:       resno = FileParserUtils::toInt(std::string(ptr + 22, 4));
        // RDKit❗✔️:     } catch (boost::bad_lexical_cast &) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Problem with residue number for PDB atom #" << serialno;
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setResidueNumber(resno);
        // RDKit❗✔️:
        // RDKit❗✔️:   double occup = 1.0;
        // RDKit❗✔️:   if (len >= 60) {
        // RDKit❗✔️:     try {
        // RDKit❗✔️:       occup = FileParserUtils::toDouble(std::string(ptr + 54, 6));
        // RDKit❗✔️:     } catch (boost::bad_lexical_cast &) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Problem with occupancy for PDB atom #" << serialno;
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setOccupancy(occup);
        // RDKit❗✔️:
        // RDKit❗✔️:   double bfactor = 0.0;
        // RDKit❗✔️:   if (len >= 66) {
        // RDKit❗✔️:     try {
        // RDKit❗✔️:       bfactor = FileParserUtils::toDouble(std::string(ptr + 60, 6));
        // RDKit❗✔️:     } catch (boost::bad_lexical_cast &) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Problem with temperature factor for PDB atom #" << serialno;
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   info->setTempFactor(bfactor);
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: void PDBBondLine(RWMol *mol, const char *ptr, unsigned int len,
        // RDKit❗✔️:                  std::map<int, Atom *> &amap, std::map<Bond *, int> &bmap) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(ptr, "bad char ptr");
        let info = AtomPdbResidueInfo::new(
            std::str::from_utf8(name).expect("validated ASCII atom name"),
            serial,
            residue.name().as_str(),
            seq.map_or(1, |s| s.seq_num()),
            chain
                .source()
                .auth_chain_id()
                .map_or_else(|| " ".to_owned(), |c| c.as_str().to_owned()),
            residue.het_flag() == Some(b'H'),
        )
        .with_alt_loc(
            row.altloc()
                .map_or_else(|| " ".to_owned(), |x| char::from(x.value()).to_string()),
        )
        .with_insertion_code(
            seq.and_then(|s| s.ins_code())
                .map_or_else(|| " ".to_owned(), |x| char::from(x).to_string()),
        )
        .with_occupancy(row.occupancy())
        .with_temp_factor(row.b_iso());
        let mut spec = AtomSpec::new(row.element())
            .with_formal_charge(row.formal_charge())
            .with_pdb_residue_info(info);
        if let Some(isotope) = row.isotope_mass_number() {
            spec = spec.with_isotope(isotope);
        }
        atoms.push(Atom::from_spec(id, spec));
        positions.push(data.coordinates().positions()[index]);
        serials.insert(serial, id);
    }
    let mut bonds: Vec<Bond> = Vec::new();
    let mut bond_map = HashMap::new();
    let mut seen = HashMap::new();
    for (&source, targets) in &data.source_state().conect_map {
        for &target in targets {
            if !crate::pdb::apply_conect_target(
                &atoms,
                &serials,
                &mut bonds,
                &mut bond_map,
                &mut seen,
                source,
                target,
            ) {
                break;
            }
        }
    }
    let mut topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())?;
    let mut coordinates = CoordinateBlock::default();
    if !positions.is_empty() {
        let is_3d = positions.iter().any(|p| p[2] != 0.0);
        coordinates
            .conformers_3d
            .push(Conformer3D::new(0, positions, is_3d));
        coordinates.source_coordinate_dim = Some(cosmolkit_model::CoordinateDimension::ThreeD);
    }
    let mut properties = MoleculeProperties::default();
    // RDKit❗✔️:   if (proximityBonding || flavor & 8) {
    // RDKit❗✔️:     StandardPDBResidueBondOrders(mol.get());
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   BasicPDBCleanup(*mol);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (sanitize) {
    // RDKit❗✔️:     if (removeHs) {
    // RDKit❗✔️:       MolOps::removeHs(*mol);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       MolOps::sanitizeMol(*mol);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // we need some properties for the chiral setup
    // RDKit❗✔️:     mol->updatePropertyCache(false);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   /* Set tetrahedral chirality from 3D co-ordinates */
    // RDKit❗✔️:   MolOps::assignChiralTypesFrom3D(*mol);
    // RDKit❗✔️:   StandardPDBResidueChirality(mol.get());
    // RDKit❗✔️:
    // RDKit❗✔️:   return mol;
    // RDKit❗✔️: }
    // RDKit❗✔️: }  // namespace
    // RDKit❗✔️:
    // RDKit❗✔️: namespace v2 {
    // RDKit❗✔️: namespace FileParsers {
    // RDKit❗✔️:
    // RDKit❗✔️: std::unique_ptr<RWMol> MolFromPDBBlock(const std::string &str,
    crate::postprocess_pdb_detached(
        &mut topology,
        &coordinates,
        crate::PdbPostprocessParams {
            proximity_bonding: params.proximity_bonding,
            flavor: params.flavor,
        },
    )?;
    if params.sanitize {
        if params.remove_hs {
            let removed = cosmolkit_core::remove_hydrogens_with_params(
                topology,
                coordinates,
                properties,
                &cosmolkit_core::RemoveHsParams::default(),
            )?;
            topology = removed.topology;
            coordinates = removed.coordinates;
            properties = removed.properties;
        } else {
            topology = cosmolkit_core::sanitize_topology(
                &topology,
                &cosmolkit_core::SanitizeParams::default(),
            )?
            .topology;
        }
    }
    let valence = cosmolkit_core::assign_valence(
        &topology,
        &cosmolkit_core::ValenceParams {
            model: cosmolkit_core::ValenceModel::RdkitLike,
            strict: false,
        },
    )?;
    let assignment = cosmolkit_core::assign_chiral_tags_from_structure(
        &topology,
        &coordinates,
        &valence,
        &cosmolkit_core::StructureTagParams::default(),
    )?;
    topology = assignment.topology;
    crate::apply_standard_pdb_residue_chirality_detached(&mut topology)?;
    Ok((topology, coordinates, properties))
}
