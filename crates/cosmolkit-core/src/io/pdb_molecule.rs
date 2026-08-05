//! RDKit-compatible PDB molecule conversion layered over `BioStructure`.

use std::collections::HashMap;

use thiserror::Error;

use crate::{
    AtomId, AtomPdbResidueInfo, AtomSpec, BioStructure, BondOrder, BondSpec, ChiralTag,
    Conformer3D, Element, Molecule, MoleculeBuildError, MoleculeBuilder, OperationError,
};

const CTD_IGNORE_H_H_CONTACTS: u32 = 0x1;
const CTD_QUICKREMOVE_H_H_CONTACTS: u32 = 0x2;
const EXTDIST: f64 = 0.45;
const MAXDIST: f64 = 5.45;
const MINDIST2: f64 = 0.16;
const MAXDIST2: f64 = 29.7025;
const HASHSIZE: usize = 1024;
const HASHMASK: i32 = 1023;
const HASHX: i32 = 571;
const HASHY: i32 = 127;
const HASHZ: i32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkitPdbMolProfile {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub flavor: u32,
    pub proximity_bonding: bool,
}

impl Default for RdkitPdbMolProfile {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hs: true,
            flavor: 0,
            proximity_bonding: true,
        }
    }
}

#[derive(Debug, Error)]
pub enum PdbMoleculeConversionError {
    #[error("BioStructure PDB read error: {0}")]
    BioRead(#[from] crate::BioReadError),
    #[error("molecule build error: {0}")]
    Build(#[from] MoleculeBuildError),
    #[error("operation error: {0}")]
    Operation(#[from] OperationError),
    #[error("stereochemistry error: {0}")]
    Stereo(#[from] crate::StereoError),
    #[error("unsupported RDKit PDB molecule conversion branch: {0}")]
    Unsupported(&'static str),
}

pub fn molecule_from_pdb_block_with_options(
    text: &str,
    profile: RdkitPdbMolProfile,
) -> Result<Molecule, PdbMoleculeConversionError> {
    let structure = BioStructure::from_pdb_str(text)?;
    bio_structure_to_rdkit_pdb_molecule(&structure, profile)
}

pub fn molecule_from_mmcif_block_with_options(
    text: &str,
    profile: RdkitPdbMolProfile,
) -> Result<Molecule, PdbMoleculeConversionError> {
    let structure = crate::io::bio::read_mmcif_atom_site_subset_from_str(text)?;
    bio_structure_to_rdkit_pdb_molecule(&structure, profile)
}

impl Molecule {
    pub fn from_pdb_block_with_options(
        text: &str,
        profile: RdkitPdbMolProfile,
    ) -> Result<Self, PdbMoleculeConversionError> {
        molecule_from_pdb_block_with_options(text, profile)
    }

    pub fn from_pdb_block(text: &str) -> Result<Self, PdbMoleculeConversionError> {
        molecule_from_pdb_block_with_options(text, RdkitPdbMolProfile::default())
    }

    pub fn from_mmcif_block_with_options(
        text: &str,
        profile: RdkitPdbMolProfile,
    ) -> Result<Self, PdbMoleculeConversionError> {
        molecule_from_mmcif_block_with_options(text, profile)
    }

    pub fn from_mmcif_block(text: &str) -> Result<Self, PdbMoleculeConversionError> {
        molecule_from_mmcif_block_with_options(text, RdkitPdbMolProfile::default())
    }
}

impl BioStructure {
    pub fn to_rdkit_pdb_molecule(
        &self,
        profile: RdkitPdbMolProfile,
    ) -> Result<Molecule, PdbMoleculeConversionError> {
        bio_structure_to_rdkit_pdb_molecule(self, profile)
    }
}

pub fn bio_structure_to_rdkit_pdb_molecule(
    structure: &BioStructure,
    profile: RdkitPdbMolProfile,
) -> Result<Molecule, PdbMoleculeConversionError> {
    let mut builder = MoleculeBuilder::new();
    let mut serial_to_atom = HashMap::<i32, AtomId>::new();
    let mut bio_atom_to_mol_atom = vec![None; structure.atoms().len()];
    let mut coords = Vec::<[f64; 3]>::new();
    let mut is_3d = false;

    for (bio_atom_index, atom) in structure.atoms().iter().enumerate() {
        let residue = &structure.residues()[atom.residue_id.index() as usize];
        if !include_atom_like_rdkit(atom, residue, structure, bio_atom_index, profile)? {
            continue;
        }
        let mol_atom = builder.add_atom(atom_spec_from_bio_atom_like_rdkit(
            atom, residue, structure,
        )?);
        if let Some(serial) = atom.source.serial {
            serial_to_atom.insert(serial.0, mol_atom);
        }
        bio_atom_to_mol_atom[bio_atom_index] = Some(mol_atom);
        let position = structure.coordinates().positions()[bio_atom_index];
        if position[2] != 0.0 {
            is_3d = true;
        }
        coords.push([
            f64::from(position[0]),
            f64::from(position[1]),
            f64::from(position[2]),
        ]);
    }

    if !coords.is_empty() {
        builder.add_conformer(Conformer3D::new(0, coords, is_3d))?;
    }

    apply_conect_records_like_rdkit(&mut builder, structure, &serial_to_atom)?;

    if profile.proximity_bonding {
        connect_the_dots_like_rdkit(&mut builder, CTD_IGNORE_H_H_CONTACTS)?;
    }
    if profile.proximity_bonding || (profile.flavor & 8) != 0 {
        standard_pdb_residue_bond_orders_like_rdkit(&mut builder)?;
    }

    basic_pdb_cleanup_like_rdkit(&mut builder);
    let topology_trust = if profile.proximity_bonding || !builder.bonds().is_empty() {
        crate::TopologyTrust::TrustedGraph
    } else {
        crate::TopologyTrust::CoordinateOnly
    };
    builder.set_topology_trust(topology_trust);
    let mut molecule = builder.build()?;
    if profile.sanitize {
        molecule = molecule.sanitize()?;
    }
    if profile.remove_hs {
        molecule = molecule.without_hydrogens()?;
    }
    // BEGIN RDKIT CPP FUNCTION parsePdbBlock chirality tail
    // RDKit✔️✔️:   /* Set tetrahedral chirality from 3D co-ordinates */
    // RDKit✔️✔️:   MolOps::assignChiralTypesFrom3D(*mol);
    // RDKit✔️✔️:   StandardPDBResidueChirality(mol.get());
    // END RDKIT CPP FUNCTION parsePdbBlock chirality tail
    crate::chemistry::stereo::assign_chiral_types_from_3d_molecule(&mut molecule, -1, true)?;
    standard_pdb_residue_chirality_like_rdkit(&mut molecule);
    Ok(molecule)
}

fn include_atom_like_rdkit(
    atom: &crate::AtomRow,
    residue: &crate::ResidueRow,
    structure: &BioStructure,
    bio_atom_index: usize,
    profile: RdkitPdbMolProfile,
) -> Result<bool, PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION PDBAtomLine filtering
    // RDKit✔️✔️:   if ((flavor & 1) == 0) {
    // RDKit✔️✔️:     // Ignore alternate locations of atoms.
    // RDKit✔️✔️:     if (len >= 17 && ptr[16] != ' ' && ptr[16] != 'A' && ptr[16] != '1') {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit❌❌:     // Ignore XPLOR pseudo atoms
    // RDKit❌❌:     if (len >= 54 && !memcmp(ptr + 30, "9999.0009999.0009999.000", 24)) {
    // RDKit❌❌:       return;
    // RDKit❌❌:     }
    // RDKit❌❌:     // Ignore NMR pseudo atoms
    // RDKit❌❌:     if (ptr[12] == ' ' && ptr[13] == 'Q') {
    // RDKit❌❌:       return;
    // RDKit❌❌:     }
    // RDKit✔️✔️:     // Ignore PDB dummy residues
    // RDKit✔️✔️:     if (len >= 20 && !memcmp(ptr + 18, "DUM", 3)) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION PDBAtomLine filtering
    if (profile.flavor & 1) != 0 {
        return Ok(true);
    }
    if atom
        .altloc
        .is_some_and(|altloc| !matches!(altloc.0, b'A' | b'1'))
    {
        return Ok(false);
    }
    if residue.name.as_str() == "DUM" {
        return Ok(false);
    }
    let position = structure.coordinates().positions()[bio_atom_index];
    if position == [9999.0, 9999.0, 9999.0] {
        return Err(PdbMoleculeConversionError::Unsupported(
            "XPLOR pseudo-atom filtering requires original PDB coordinate field text",
        ));
    }
    Ok(true)
}

fn atom_spec_from_bio_atom_like_rdkit(
    atom: &crate::AtomRow,
    residue: &crate::ResidueRow,
    structure: &BioStructure,
) -> Result<AtomSpec, PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION PDBAtomLine atom state
    // RDKit❗✔️:   Atom *atom = (Atom *)nullptr;
    // RDKit❗✔️:   char symb[3];
    // RDKit❌❌:   // Attempt #1:  Atomic Symbol in columns 76 and 77
    // RDKit❌❌:   // Attempt #2: Atomic Symbol from PDB atom name
    // RDKit✔️✔️:   if (charge != 0) {
    // RDKit✔️✔️:     atom->setFormalCharge(charge);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   tmp = std::string(ptr + 12, 4);
    // RDKit✔️✔️:   AtomPDBResidueInfo *info = new AtomPDBResidueInfo(tmp, serialno);
    // RDKit✔️✔️:   atom->setMonomerInfo(info);
    // RDKit✔️✔️:   info->setResidueName(tmp);
    // RDKit✔️✔️:   if (ptr[0] == 'H') {
    // RDKit✔️✔️:     info->setIsHeteroAtom(true);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   info->setAltLoc(tmp);
    // RDKit✔️✔️:   info->setChainId(tmp);
    // RDKit✔️✔️:   info->setInsertionCode(tmp);
    // RDKit✔️✔️:   info->setResidueNumber(resno);
    // RDKit✔️✔️:   info->setOccupancy(occup);
    // RDKit✔️✔️:   info->setTempFactor(bfactor);
    // END RDKIT CPP FUNCTION PDBAtomLine atom state
    let mut spec = AtomSpec::new(atom.element);
    if let Some(charge) = atom.formal_charge {
        spec = spec.with_formal_charge(charge);
    }
    if atom.element == Element::H {
        let atom_name = atom_name_string(atom.name);
        if atom_name.trim() == "D" {
            spec = spec.with_isotope(2);
        } else if atom_name.trim() == "T" {
            spec = spec.with_isotope(3);
        }
    }
    let serial = atom.source.serial.map_or(0, |serial| serial.0);
    let chain = &structure.chains()[residue.chain_id.index() as usize];
    let chain_id = chain
        .source
        .auth_chain_id
        .map_or_else(|| " ".to_string(), |chain_id| chain_id.as_str().to_string());
    let residue_number = residue.source.seq_id.map_or(1, |seq_id| seq_id.seq_num);
    let insertion_code = residue
        .source
        .seq_id
        .and_then(|seq_id| seq_id.ins_code)
        .map_or_else(|| " ".to_string(), |code| char::from(code).to_string());
    let alt_loc = atom.altloc.map_or_else(
        || " ".to_string(),
        |altloc| char::from(altloc.0).to_string(),
    );
    let is_hetero_atom = residue.het_flag == Some('H');
    let pdb_info = AtomPdbResidueInfo::new(
        atom_name_string(atom.name),
        serial,
        residue.name.as_str().to_string(),
        residue_number,
        chain_id,
        is_hetero_atom,
    )
    .with_alt_loc(alt_loc)
    .with_insertion_code(insertion_code)
    .with_occupancy(f64::from(atom.occupancy.unwrap_or(1.0)))
    .with_temp_factor(f64::from(atom.b_iso.unwrap_or(0.0)));
    Ok(spec.with_pdb_residue_info(pdb_info))
}

fn apply_conect_records_like_rdkit(
    builder: &mut MoleculeBuilder,
    structure: &BioStructure,
    serial_to_atom: &HashMap<i32, AtomId>,
) -> Result<(), PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION PDBBondLine
    // RDKit✔️✔️:   if (len < 16) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string tmp(ptr + 6, 5);
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     src = FileParserUtils::toInt(tmp);
    // RDKit✔️✔️:     if (amap.find(src) == amap.end()) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int pos = 11; pos + 5 <= len; pos += 5) {
    // RDKit✔️✔️:     if (!memcmp(ptr + pos, "     ", 5)) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     dst = FileParserUtils::toInt(std::string(ptr + pos, 5));
    // RDKit✔️✔️:     if (dst == src || amap.find(dst) == amap.end()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // END RDKIT CPP FUNCTION PDBBondLine
    let mut bond_seen = HashMap::<crate::BondId, u8>::new();
    let mut sorted_sources = structure.conect_map.keys().copied().collect::<Vec<_>>();
    sorted_sources.sort_unstable();
    for src in sorted_sources {
        let Some(&begin) = serial_to_atom.get(&src) else {
            continue;
        };
        for dst in &structure.conect_map[&src] {
            if *dst == src {
                continue;
            }
            let Some(&end) = serial_to_atom.get(dst) else {
                continue;
            };
            apply_one_conect_bond_like_rdkit(builder, begin, end, src, *dst, &mut bond_seen)?;
        }
    }
    Ok(())
}

fn apply_one_conect_bond_like_rdkit(
    builder: &mut MoleculeBuilder,
    begin: AtomId,
    end: AtomId,
    src_serial: i32,
    dst_serial: i32,
    bond_seen: &mut HashMap<crate::BondId, u8>,
) -> Result<(), PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION PDBBondLine bond order bitmap
    // RDKit✔️✔️:       Bond *bond =
    // RDKit✔️✔️:           mol->getBondBetweenAtoms(amap[src]->getIdx(), amap[dst]->getIdx());
    // RDKit✔️✔️:       if (bond && bond->getBondType() != Bond::ZERO) {
    // RDKit✔️✔️:         // Here we use a single byte bitmap to count duplicates
    // RDKit✔️✔️:         // Low nibble counts src < dst, high nibble for src > dst
    // RDKit✔️✔️:         int seen = bmap[bond];
    // RDKit✔️✔️:         if (src < dst) {
    // RDKit✔️✔️:           if ((seen & 0x0f) == 0x01) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x02;
    // RDKit✔️✔️:             if ((seen & 0x20) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if ((seen & 0x0f) == 0x03) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x04;
    // RDKit✔️✔️:             if ((seen & 0x40) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::TRIPLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if ((seen & 0x0f) == 0x07) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x08;
    // RDKit✔️✔️:             if ((seen & 0x80) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::QUADRUPLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         } else /* src < dst */ {
    // RDKit✔️✔️:           if ((seen & 0xf0) == 0x10) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x20;
    // RDKit✔️✔️:             if ((seen & 0x02) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if ((seen & 0xf0) == 0x30) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x40;
    // RDKit✔️✔️:             if ((seen & 0x04) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::TRIPLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if ((seen & 0xf0) == 0x70) {
    // RDKit✔️✔️:             bmap[bond] = seen | 0x80;
    // RDKit✔️✔️:             if ((seen & 0x08) == 0) {
    // RDKit✔️✔️:               bond->setBondType(Bond::QUADRUPLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (!bond) {
    // RDKit❌❌:         if (IsBlacklistedPair(amap[src], amap[dst])) {
    // RDKit❌❌:           bond = new Bond(Bond::ZERO);
    // RDKit❌❌:         } else {
    // RDKit✔️✔️:           bond = new Bond(Bond::SINGLE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         mol->addBond(bond, true);
    // RDKit✔️✔️:         bmap[bond] = (src < dst) ? 0x01 : 0x10;
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION PDBBondLine bond order bitmap
    let Some(bond) = builder.bond_between_atoms(begin, end) else {
        let order = if is_blacklisted_pair_like_rdkit(builder, begin, end) {
            BondOrder::Zero
        } else {
            BondOrder::Single
        };
        let bond = builder.add_bond(BondSpec::new(begin, end, order))?;
        bond_seen.insert(bond, if src_serial < dst_serial { 0x01 } else { 0x10 });
        return Ok(());
    };
    if builder
        .bond(bond)
        .is_some_and(|bond| bond.order() == BondOrder::Zero)
    {
        return Ok(());
    }
    let seen = *bond_seen.get(&bond).unwrap_or(&0);
    let mut next_seen = seen;
    let mut next_order = None;
    if src_serial < dst_serial {
        if (seen & 0x0f) == 0x01 {
            next_seen = seen | 0x02;
            if (seen & 0x20) == 0 {
                next_order = Some(BondOrder::Double);
            }
        } else if (seen & 0x0f) == 0x03 {
            next_seen = seen | 0x04;
            if (seen & 0x40) == 0 {
                next_order = Some(BondOrder::Triple);
            }
        } else if (seen & 0x0f) == 0x07 {
            next_seen = seen | 0x08;
            if (seen & 0x80) == 0 {
                next_order = Some(BondOrder::Quadruple);
            }
        }
    } else if (seen & 0xf0) == 0x10 {
        next_seen = seen | 0x20;
        if (seen & 0x02) == 0 {
            next_order = Some(BondOrder::Double);
        }
    } else if (seen & 0xf0) == 0x30 {
        next_seen = seen | 0x40;
        if (seen & 0x04) == 0 {
            next_order = Some(BondOrder::Triple);
        }
    } else if (seen & 0xf0) == 0x70 {
        next_seen = seen | 0x80;
        if (seen & 0x08) == 0 {
            next_order = Some(BondOrder::Quadruple);
        }
    }
    bond_seen.insert(bond, next_seen);
    if let Some(order) = next_order {
        builder.set_bond_order(bond, order)?;
    }
    Ok(())
}

fn same_pdb_residue_like_rdkit(left: &AtomPdbResidueInfo, right: &AtomPdbResidueInfo) -> bool {
    // BEGIN RDKIT CPP FUNCTION SamePDBResidue
    // RDKit✔️✔️: bool SamePDBResidue(AtomPDBResidueInfo *p, AtomPDBResidueInfo *q) {
    // RDKit✔️✔️:   return p->getResidueNumber() == q->getResidueNumber() &&
    // RDKit✔️✔️:          p->getResidueName() == q->getResidueName() &&
    // RDKit✔️✔️:          p->getChainId() == q->getChainId() &&
    // RDKit✔️✔️:          p->getInsertionCode() == q->getInsertionCode();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SamePDBResidue
    left.residue_number() == right.residue_number()
        && left.residue_name() == right.residue_name()
        && left.chain_id() == right.chain_id()
        && left.insertion_code() == right.insertion_code()
}

fn is_blacklisted_atom_like_rdkit(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION IsBlacklistedAtom
    // RDKit✔️✔️: static bool IsBlacklistedAtom(Atom *atom) {
    // RDKit✔️✔️:   // blacklist metals, noble gasses and halogens
    // RDKit✔️✔️:   int elem = atom->getAtomicNum();
    // RDKit✔️✔️:   // make an inverse query (non-metals and metaloids)
    // RDKit✔️✔️:   return !((5 <= elem && elem <= 8) || (14 <= elem && elem <= 16) ||
    // RDKit✔️✔️:            (32 <= elem && elem <= 34) || (51 <= elem && elem <= 52));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION IsBlacklistedAtom
    !((5..=8).contains(&atomic_number)
        || (14..=16).contains(&atomic_number)
        || (32..=34).contains(&atomic_number)
        || (51..=52).contains(&atomic_number))
}

fn is_blacklisted_pair_like_rdkit(builder: &MoleculeBuilder, begin: AtomId, end: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION IsBlacklistedPair
    // RDKit✔️✔️: bool IsBlacklistedPair(Atom *beg_atom, Atom *end_atom) {
    // RDKit✔️✔️:   PRECONDITION(beg_atom, "empty atom");
    // RDKit✔️✔️:   PRECONDITION(end_atom, "empty atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *beg_info = (AtomPDBResidueInfo *)beg_atom->getMonomerInfo();
    // RDKit✔️✔️:   auto *end_info = (AtomPDBResidueInfo *)end_atom->getMonomerInfo();
    // RDKit✔️✔️:   if (!beg_info || beg_info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!end_info || end_info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!SamePDBResidue(beg_info, end_info)) {
    // RDKit✔️✔️:     if (IsBlacklistedAtom(beg_atom) || IsBlacklistedAtom(end_atom)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Dont make bonds to waters
    // RDKit✔️✔️:     if (beg_info->getResidueName() == "HOH" ||
    // RDKit✔️✔️:         end_info->getResidueName() == "HOH") {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION IsBlacklistedPair
    let Some(begin_atom) = builder.atoms().get(begin.index()) else {
        return false;
    };
    let Some(end_atom) = builder.atoms().get(end.index()) else {
        return false;
    };
    let Some(begin_info) = begin_atom.pdb_residue_info() else {
        return false;
    };
    let Some(end_info) = end_atom.pdb_residue_info() else {
        return false;
    };
    if !same_pdb_residue_like_rdkit(begin_info, end_info) {
        if is_blacklisted_atom_like_rdkit(begin_atom.atomic_number())
            || is_blacklisted_atom_like_rdkit(end_atom.atomic_number())
        {
            return true;
        }
        if begin_info.residue_name() == "HOH" || end_info.residue_name() == "HOH" {
            return true;
        }
    }
    false
}

#[derive(Debug, Clone, Copy)]
struct ProximityEntry {
    x: f32,
    y: f32,
    z: f32,
    r: f32,
    atm: usize,
    hash: i32,
    next: i32,
    elem: u8,
}

fn rdkit_covalent_radius(atomic_number: u8) -> Result<f32, PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRcovalent / atomicData::Rcov
    // RDKit✔️✔️: double getRcovalent(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Rcov();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double Rcov() const { return rCov; }
    // RDKit✔️✔️: //  rCov (https://doi.org/10.1039/B801115J). 1.9 if unknown. the low spin value
    // RDKit✔️✔️: //  was taken for transition metals
    // RDKit✔️✔️: const std::string periodicTableAtomData =
    // RDKit✔️✔️:     R"DAT(0	*	0	0	0	0	0	0	0	0	-1
    // RDKit✔️✔️: 1	H		1	0.31	0.33	1.2	1.008	1	1	1.007825032	1
    // RDKit✔️✔️: 2	He	1	0.28	0.7	1.4	4.003	2	4	4.002603254	0
    // RDKit✔️✔️: 3	Li	2	1.28	1.23	2.2	6.941	1	7	7.01600455	1	-1
    // RDKit✔️✔️: 4	Be	2	0.96	0.9	1.9	9.012	2	9	9.0121822	2
    // RDKit✔️✔️: 5	B	2	0.84	0.82	1.8	10.812	3	11	11.0093054	3
    // RDKit✔️✔️: 6	C	2	0.76	0.77	1.7	12.011	4	12	12	4
    // RDKit✔️✔️: 7	N	2	0.71	0.7	1.6	14.007	5	14	14.003074	3
    // RDKit✔️✔️: 8	O	2	0.66	0.66	1.55	15.999	6	16	15.99491462	2
    // RDKit✔️✔️: 9	F	2	0.57	0.611	1.5	18.998	7	19	18.99840322	1
    // RDKit✔️✔️: 10	Ne	2	0.58	0.7	1.54	20.18	8	20	19.99244018	0
    // RDKit❗✔️: // Full rCov column through Og is stored in RDKIT_COVALENT_RADII below.
    // END RDKIT CPP FUNCTION PeriodicTable::getRcovalent / atomicData::Rcov
    RDKIT_COVALENT_RADII
        .get(usize::from(atomic_number))
        .copied()
        .ok_or(PdbMoleculeConversionError::Unsupported(
            "RDKit covalent radius lookup requires an atomic number in periodicTableAtomData",
        ))
}

const RDKIT_COVALENT_RADII: [f32; 119] = [
    0.0, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07,
    1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.52, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22,
    1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45,
    1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98,
    1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
    1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.6, 2.2, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80,
    1.69, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.36,
    1.43, 1.62, 1.75, 1.65, 1.57,
];

fn is_bonded_like_rdkit(p: &ProximityEntry, q: &ProximityEntry, flags: u32) -> bool {
    // BEGIN RDKIT CPP FUNCTION IsBonded
    // RDKit✔️✔️: static bool IsBonded(ProximityEntry *p, ProximityEntry *q, unsigned int flags) {
    // RDKit✔️✔️:   if (flags & ctdIGNORE_H_H_CONTACTS && p->elem == 1 && q->elem == 1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double dx = (double)p->x - (double)q->x;
    // RDKit✔️✔️:   double dist2 = dx * dx;
    // RDKit✔️✔️:   if (dist2 > MAXDIST2) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double dy = (double)p->y - (double)q->y;
    // RDKit✔️✔️:   dist2 += dy * dy;
    // RDKit✔️✔️:   if (dist2 > MAXDIST2) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double dz = (double)p->z - (double)q->z;
    // RDKit✔️✔️:   dist2 += dz * dz;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (dist2 > MAXDIST2 || dist2 < MINDIST2) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double radius = (double)p->r + (double)q->r + EXTDIST;
    // RDKit✔️✔️:   return dist2 <= radius * radius;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION IsBonded
    if flags & CTD_IGNORE_H_H_CONTACTS != 0 && p.elem == 1 && q.elem == 1 {
        return false;
    }
    let dx = f64::from(p.x) - f64::from(q.x);
    let mut dist2 = dx * dx;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dy = f64::from(p.y) - f64::from(q.y);
    dist2 += dy * dy;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dz = f64::from(p.z) - f64::from(q.z);
    dist2 += dz * dz;
    if !(MINDIST2..=MAXDIST2).contains(&dist2) {
        return false;
    }
    let radius = f64::from(p.r) + f64::from(q.r) + EXTDIST;
    dist2 <= radius * radius
}

fn other_atom_for_bond(bond: &crate::Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn distance3(a: [f64; 3], b: [f64; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt() as f32
}

fn cleanup_multivalent_hydrogens_like_rdkit(builder: &mut MoleculeBuilder, flags: u32) {
    // BEGIN RDKIT CPP FUNCTION ConnectTheDots_Large cleanup pass
    // RDKit✔️✔️:   // Cleanup pass
    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     Atom *atom = mol->getAtomWithIdx(i);
    // RDKit✔️✔️:     unsigned int elem = atom->getAtomicNum();
    // RDKit✔️✔️:     // detect multivalent Hs, which could happen with ConnectTheDots
    // RDKit✔️✔️:     if (elem == 1 && atom->getDegree() > 1) {
    // RDKit✔️✔️:       // if there's an H neighbor and a non-H neighbor, remove the bond to the H
    // RDKit✔️✔️:       if (flags & ctdQUICKREMOVE_H_H_CONTACTS) {
    // RDKit✔️✔️:         Bond *bondToH = nullptr;
    // RDKit✔️✔️:         Bond *bondToNonH = nullptr;
    // RDKit✔️✔️:         for (auto bond : mol->atomBonds(atom)) {
    // RDKit✔️✔️:           if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
    // RDKit✔️✔️:             bondToH = bond;
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             bondToNonH = bond;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (bondToH && bondToNonH) {
    // RDKit✔️✔️:           mol->removeBond(bondToH->getBeginAtomIdx(), bondToH->getEndAtomIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // if we now have a degree of 1, we're done
    // RDKit✔️✔️:       if (atom->getDegree() == 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       auto *atom_info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
    // RDKit✔️✔️:       // cut all but shortest Bond
    // RDKit✔️✔️:       RDGeom::Point3D p = conf->getAtomPos(i);
    // RDKit✔️✔️:       RDKit::RWMol::ADJ_ITER nbr, end_nbr;
    // RDKit✔️✔️:       boost::tie(nbr, end_nbr) = mol->getAtomNeighbors(atom);
    // RDKit✔️✔️:       float best = 10000;
    // RDKit✔️✔️:       unsigned int best_idx = mol->getNumAtoms() + 1;
    // RDKit✔️✔️:       while (nbr != end_nbr) {
    // RDKit✔️✔️:         RDGeom::Point3D pn = conf->getAtomPos(*nbr);
    // RDKit✔️✔️:         float d = (p - pn).length();
    // RDKit✔️✔️:         auto *n_info =
    // RDKit✔️✔️:             (AtomPDBResidueInfo *)(mol->getAtomWithIdx(*nbr)->getMonomerInfo());
    // RDKit✔️✔️:         if (d < best &&
    // RDKit✔️✔️:             ((!atom_info || !n_info) ||
    // RDKit✔️✔️:              atom_info->getResidueNumber() == n_info->getResidueNumber())) {
    // RDKit✔️✔️:           best = d;
    // RDKit✔️✔️:           best_idx = *nbr;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++nbr;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // iterate again and remove all but closest
    // RDKit✔️✔️:       boost::tie(nbr, end_nbr) = mol->getAtomNeighbors(atom);
    // RDKit✔️✔️:       while (nbr != end_nbr) {
    // RDKit✔️✔️:         if (*nbr == best_idx) {
    // RDKit✔️✔️:           Bond *bond = mol->getBondBetweenAtoms(i, *nbr);
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);  // make sure this one is single
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           mol->removeBond(i, *nbr);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++nbr;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION ConnectTheDots_Large cleanup pass
    let Some(conformer) = builder.conformers_3d().first() else {
        return;
    };
    let coords = conformer.coordinates().to_vec();
    let count = builder.atoms().len();
    for i in 0..count {
        let atom_id = AtomId::new(i);
        if builder.atoms()[i].atomic_number() != 1 || builder.degree(atom_id) <= 1 {
            continue;
        }
        if flags & CTD_QUICKREMOVE_H_H_CONTACTS != 0 {
            let mut bond_to_h = None;
            let mut bond_to_non_h = None;
            for bond_id in builder.neighbor_bonds(atom_id) {
                let Some(bond) = builder.bond(*bond_id) else {
                    continue;
                };
                let other = other_atom_for_bond(bond, atom_id);
                if builder.atoms()[other.index()].atomic_number() == 1 {
                    bond_to_h = Some(*bond_id);
                } else {
                    bond_to_non_h = Some(*bond_id);
                }
            }
            if bond_to_h.is_some()
                && bond_to_non_h.is_some()
                && let Some(bond) = bond_to_h.and_then(|bond_id| builder.bond(bond_id).cloned())
            {
                builder.remove_bond_between_atoms(bond.begin(), bond.end());
            }
        }
        if builder.degree(atom_id) == 1 {
            continue;
        }

        let atom_info = builder.atoms()[i].pdb_residue_info();
        let p = coords[i];
        let mut best = 10000.0_f32;
        let mut best_idx = count + 1;
        let neighbors = builder.neighbor_bonds(atom_id).to_vec();
        for bond_id in &neighbors {
            let Some(bond) = builder.bond(*bond_id) else {
                continue;
            };
            let other = other_atom_for_bond(bond, atom_id);
            let d = distance3(p, coords[other.index()]);
            let n_info = builder.atoms()[other.index()].pdb_residue_info();
            if d < best
                && (atom_info.is_none()
                    || n_info.is_none()
                    || atom_info.unwrap().residue_number() == n_info.unwrap().residue_number())
            {
                best = d;
                best_idx = other.index();
            }
        }
        let neighbors = builder.neighbor_bonds(atom_id).to_vec();
        for bond_id in neighbors {
            let Some(bond) = builder.bond(bond_id).cloned() else {
                continue;
            };
            let other = other_atom_for_bond(&bond, atom_id);
            if other.index() == best_idx {
                let _ = builder.set_bond_order(bond_id, BondOrder::Single);
            } else {
                builder.remove_bond_between_atoms(atom_id, other);
            }
        }
    }
}

fn connect_the_dots_like_rdkit(
    builder: &mut MoleculeBuilder,
    flags: u32,
) -> Result<(), PdbMoleculeConversionError> {
    // BEGIN RDKIT CPP FUNCTION ConnectTheDots / ConnectTheDots_Large
    // RDKit✔️✔️: void ConnectTheDots(RWMol *mol, unsigned int flags) {
    // RDKit✔️✔️:   if (!mol || !mol->getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // Determine optimal algorithm to use by getNumAtoms()?
    // RDKit✔️✔️:   ConnectTheDots_Large(mol, flags);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: static void ConnectTheDots_Large(RWMol *mol, unsigned int flags) {
    // RDKit✔️✔️:   int HashTable[HASHSIZE];
    // RDKit✔️✔️:   memset(HashTable, -1, sizeof(HashTable));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int count = mol->getNumAtoms();
    // RDKit✔️✔️:   auto *tmp = (ProximityEntry *)malloc(count * sizeof(ProximityEntry));
    // RDKit✔️✔️:   CHECK_INVARIANT(tmp, "bad allocation");
    // RDKit✔️✔️:   PeriodicTable *table = PeriodicTable::getTable();
    // RDKit✔️✔️:   Conformer *conf = &mol->getConformer();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     Atom *atom = mol->getAtomWithIdx(i);
    // RDKit✔️✔️:     unsigned int elem = atom->getAtomicNum();
    // RDKit✔️✔️:     RDGeom::Point3D p = conf->getAtomPos(i);
    // RDKit✔️✔️:     ProximityEntry *tmpi = tmp + i;
    // RDKit✔️✔️:     tmpi->x = (float)p.x;
    // RDKit✔️✔️:     tmpi->y = (float)p.y;
    // RDKit✔️✔️:     tmpi->z = (float)p.z;
    // RDKit✔️✔️:     tmpi->r = (float)table->getRcovalent(elem);
    // RDKit✔️✔️:     tmpi->atm = i;
    // RDKit✔️✔️:     tmpi->elem = elem;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int hash = HASHX * (int)(p.x / MAXDIST) + HASHY * (int)(p.y / MAXDIST) +
    // RDKit✔️✔️:                HASHZ * (int)(p.z / MAXDIST);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (int dx = -HASHX; dx <= HASHX; dx += HASHX) {
    // RDKit✔️✔️:       for (int dy = -HASHY; dy <= HASHY; dy += HASHY) {
    // RDKit✔️✔️:         for (int dz = -HASHZ; dz <= HASHZ; dz += HASHZ) {
    // RDKit✔️✔️:           int probe = hash + dx + dy + dz;
    // RDKit✔️✔️:           int list = HashTable[probe & HASHMASK];
    // RDKit✔️✔️:           while (list != -1) {
    // RDKit✔️✔️:             ProximityEntry *tmpj = &tmp[list];
    // RDKit✔️✔️:             if (tmpj->hash == probe && IsBonded(tmpi, tmpj, flags) &&
    // RDKit✔️✔️:                 !mol->getBondBetweenAtoms(tmpi->atm, tmpj->atm) &&
    // RDKit✔️✔️:                 !IsBlacklistedPair(atom, mol->getAtomWithIdx(tmpj->atm))) {
    // RDKit✔️✔️:               mol->addBond(tmpi->atm, tmpj->atm, Bond::SINGLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             list = tmpj->next;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     int list = hash & HASHMASK;
    // RDKit✔️✔️:     tmpi->next = HashTable[list];
    // RDKit✔️✔️:     HashTable[list] = i;
    // RDKit✔️✔️:     tmpi->hash = hash;
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   // Cleanup pass body is reproduced in cleanup_multivalent_hydrogens_like_rdkit().
    // RDKit✔️✔️:   free(tmp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ConnectTheDots / ConnectTheDots_Large
    let Some(conformer) = builder.conformers_3d().first() else {
        return Ok(());
    };
    let coords = conformer.coordinates().to_vec();
    let mut hash_table = [-1_i32; HASHSIZE];
    let count = builder.atoms().len();
    let mut tmp = Vec::<ProximityEntry>::with_capacity(count);
    for i in 0..count {
        let atom = &builder.atoms()[i];
        let elem = atom.atomic_number();
        let p = coords[i];
        let hash = HASHX * (p[0] / MAXDIST) as i32
            + HASHY * (p[1] / MAXDIST) as i32
            + HASHZ * (p[2] / MAXDIST) as i32;
        let tmpi = ProximityEntry {
            x: p[0] as f32,
            y: p[1] as f32,
            z: p[2] as f32,
            r: rdkit_covalent_radius(elem)?,
            atm: i,
            hash,
            next: -1,
            elem,
        };
        for dx in [-HASHX, 0, HASHX] {
            for dy in [-HASHY, 0, HASHY] {
                for dz in [-HASHZ, 0, HASHZ] {
                    let probe = hash + dx + dy + dz;
                    let mut list = hash_table[(probe & HASHMASK) as usize];
                    while list != -1 {
                        let tmpj = tmp[list as usize];
                        let begin = AtomId::new(tmpi.atm);
                        let end = AtomId::new(tmpj.atm);
                        if tmpj.hash == probe
                            && is_bonded_like_rdkit(&tmpi, &tmpj, flags)
                            && builder.bond_between_atoms(begin, end).is_none()
                            && !is_blacklisted_pair_like_rdkit(builder, begin, end)
                        {
                            builder.add_bond(BondSpec::new(begin, end, BondOrder::Single))?;
                        }
                        list = tmpj.next;
                    }
                }
            }
        }
        let list = (hash & HASHMASK) as usize;
        let mut stored = tmpi;
        stored.next = hash_table[list];
        hash_table[list] = i as i32;
        tmp.push(stored);
    }
    cleanup_multivalent_hydrogens_like_rdkit(builder, flags);
    Ok(())
}

const fn bcnam(a: u8, b: u8, c: u8) -> u32 {
    ((a as u32) << 16) | ((b as u32) << 8) | c as u32
}

const fn bcatm(a: u8, b: u8, c: u8, d: u8) -> u32 {
    ((a as u32) << 24) | ((b as u32) << 16) | ((c as u32) << 8) | d as u32
}

fn fixed_width_code(bytes: &[u8], width: usize) -> u32 {
    let mut code = 0_u32;
    for index in 0..width {
        let byte = bytes.get(index).copied().unwrap_or(b'\0');
        code = (code << 8) | u32::from(byte);
    }
    code
}

fn standard_pdb_double_bond_table_like_rdkit(
    mut rescode: u32,
    mut atm1: u32,
    mut atm2: u32,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION StandardPDBDoubleBond table
    // RDKit✔️✔️: static bool StandardPDBDoubleBond(unsigned int rescode, unsigned int atm1,
    // RDKit✔️✔️:                                   unsigned int atm2) {
    // RDKit✔️✔️:   if (atm1 > atm2) {
    // RDKit✔️✔️:     unsigned int tmp = atm1;
    // RDKit✔️✔️:     atm1 = atm2;
    // RDKit✔️✔️:     atm2 = tmp;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   switch (rescode) {
    // RDKit✔️✔️:     case BCNAM('A', 'C', 'E'):
    // RDKit✔️✔️:     case BCNAM('A', 'L', 'A'):
    // RDKit✔️✔️:     case BCNAM('C', 'Y', 'S'):
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'Y'):
    // RDKit✔️✔️:     case BCNAM('I', 'L', 'E'):
    // RDKit✔️✔️:     case BCNAM('L', 'E', 'U'):
    // RDKit✔️✔️:     case BCNAM('L', 'Y', 'S'):
    // RDKit✔️✔️:     case BCNAM('M', 'E', 'T'):
    // RDKit✔️✔️:     case BCNAM('P', 'R', 'O'):
    // RDKit✔️✔️:     case BCNAM('S', 'E', 'R'):
    // RDKit✔️✔️:     case BCNAM('T', 'H', 'R'):
    // RDKit✔️✔️:     case BCNAM('V', 'A', 'L'):
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', ' ', ' ')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case BCNAM('A', 'R', 'G'):
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', ' ', ' ')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', 'Z', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'N', 'H', '2')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case BCNAM('A', 'S', 'N'):
    // RDKit✔️✔️:     case BCNAM('A', 'S', 'P'):
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', ' ', ' ')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', 'G', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', 'D', '1')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'N'):
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'U'):
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', ' ', ' ')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atm1 == BCATM(' ', 'C', 'D', ' ') &&
    // RDKit✔️✔️:           atm2 == BCATM(' ', 'O', 'E', '1')) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit❗✔️:     // Remaining source table rows are encoded in `STANDARD_PDB_DOUBLE_BONDS`.
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StandardPDBDoubleBond table
    if atm1 > atm2 {
        std::mem::swap(&mut atm1, &mut atm2);
    }
    STANDARD_PDB_DOUBLE_BONDS
        .iter()
        .any(|&(res, left, right)| res == rescode && left == atm1 && right == atm2)
}

const STANDARD_PDB_DOUBLE_BONDS: &[(u32, u32, u32)] = &[
    (
        bcnam(b'A', b'C', b'E'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'A', b'L', b'A'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'C', b'Y', b'S'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'G', b'L', b'Y'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'I', b'L', b'E'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'L', b'E', b'U'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'L', b'Y', b'S'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'M', b'E', b'T'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'P', b'R', b'O'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'S', b'E', b'R'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'T', b'H', b'R'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'V', b'A', b'L'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'A', b'R', b'G'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'A', b'R', b'G'),
        bcatm(b' ', b'C', b'Z', b' '),
        bcatm(b' ', b'N', b'H', b'2'),
    ),
    (
        bcnam(b'A', b'S', b'N'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'A', b'S', b'N'),
        bcatm(b' ', b'C', b'G', b' '),
        bcatm(b' ', b'O', b'D', b'1'),
    ),
    (
        bcnam(b'A', b'S', b'P'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'A', b'S', b'P'),
        bcatm(b' ', b'C', b'G', b' '),
        bcatm(b' ', b'O', b'D', b'1'),
    ),
    (
        bcnam(b'G', b'L', b'N'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'G', b'L', b'N'),
        bcatm(b' ', b'C', b'D', b' '),
        bcatm(b' ', b'O', b'E', b'1'),
    ),
    (
        bcnam(b'G', b'L', b'U'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'G', b'L', b'U'),
        bcatm(b' ', b'C', b'D', b' '),
        bcatm(b' ', b'O', b'E', b'1'),
    ),
    (
        bcnam(b'H', b'I', b'S'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'H', b'I', b'S'),
        bcatm(b' ', b'C', b'D', b'2'),
        bcatm(b' ', b'C', b'G', b' '),
    ),
    (
        bcnam(b'H', b'I', b'S'),
        bcatm(b' ', b'C', b'E', b'1'),
        bcatm(b' ', b'N', b'D', b'1'),
    ),
    (
        bcnam(b'P', b'H', b'E'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'P', b'H', b'E'),
        bcatm(b' ', b'C', b'D', b'1'),
        bcatm(b' ', b'C', b'G', b' '),
    ),
    (
        bcnam(b'P', b'H', b'E'),
        bcatm(b' ', b'C', b'D', b'2'),
        bcatm(b' ', b'C', b'E', b'2'),
    ),
    (
        bcnam(b'P', b'H', b'E'),
        bcatm(b' ', b'C', b'E', b'1'),
        bcatm(b' ', b'C', b'Z', b' '),
    ),
    (
        bcnam(b'T', b'Y', b'R'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'T', b'Y', b'R'),
        bcatm(b' ', b'C', b'D', b'1'),
        bcatm(b' ', b'C', b'G', b' '),
    ),
    (
        bcnam(b'T', b'Y', b'R'),
        bcatm(b' ', b'C', b'D', b'2'),
        bcatm(b' ', b'C', b'E', b'2'),
    ),
    (
        bcnam(b'T', b'Y', b'R'),
        bcatm(b' ', b'C', b'E', b'1'),
        bcatm(b' ', b'C', b'Z', b' '),
    ),
    (
        bcnam(b'T', b'R', b'P'),
        bcatm(b' ', b'C', b' ', b' '),
        bcatm(b' ', b'O', b' ', b' '),
    ),
    (
        bcnam(b'T', b'R', b'P'),
        bcatm(b' ', b'C', b'D', b'1'),
        bcatm(b' ', b'C', b'G', b' '),
    ),
    (
        bcnam(b'T', b'R', b'P'),
        bcatm(b' ', b'C', b'D', b'2'),
        bcatm(b' ', b'C', b'E', b'2'),
    ),
    (
        bcnam(b'T', b'R', b'P'),
        bcatm(b' ', b'C', b'E', b'3'),
        bcatm(b' ', b'C', b'Z', b'3'),
    ),
    (
        bcnam(b'T', b'R', b'P'),
        bcatm(b' ', b'C', b'H', b'2'),
        bcatm(b' ', b'C', b'Z', b'2'),
    ),
    (
        bcnam(b' ', b' ', b'A'),
        bcatm(b' ', b'C', b'6', b' '),
        bcatm(b' ', b'N', b'1', b' '),
    ),
    (
        bcnam(b' ', b' ', b'A'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b' ', b'A'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'C', b'5', b' '),
    ),
    (
        bcnam(b' ', b' ', b'A'),
        bcatm(b' ', b'C', b'8', b' '),
        bcatm(b' ', b'N', b'7', b' '),
    ),
    (
        bcnam(b' ', b' ', b'A'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b'D', b'A'),
        bcatm(b' ', b'C', b'6', b' '),
        bcatm(b' ', b'N', b'1', b' '),
    ),
    (
        bcnam(b' ', b'D', b'A'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b'D', b'A'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'C', b'5', b' '),
    ),
    (
        bcnam(b' ', b'D', b'A'),
        bcatm(b' ', b'C', b'8', b' '),
        bcatm(b' ', b'N', b'7', b' '),
    ),
    (
        bcnam(b' ', b'D', b'A'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b' ', b'G'),
        bcatm(b' ', b'C', b'6', b' '),
        bcatm(b' ', b'O', b'6', b' '),
    ),
    (
        bcnam(b' ', b' ', b'G'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b' ', b'G'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'C', b'5', b' '),
    ),
    (
        bcnam(b' ', b' ', b'G'),
        bcatm(b' ', b'C', b'8', b' '),
        bcatm(b' ', b'N', b'7', b' '),
    ),
    (
        bcnam(b' ', b' ', b'G'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b'D', b'G'),
        bcatm(b' ', b'C', b'6', b' '),
        bcatm(b' ', b'O', b'6', b' '),
    ),
    (
        bcnam(b' ', b'D', b'G'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b'D', b'G'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'C', b'5', b' '),
    ),
    (
        bcnam(b' ', b'D', b'G'),
        bcatm(b' ', b'C', b'8', b' '),
        bcatm(b' ', b'N', b'7', b' '),
    ),
    (
        bcnam(b' ', b'D', b'G'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b' ', b'C'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b' ', b'C'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b' ', b'C'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b' ', b'C'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b'D', b'C'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b'D', b'C'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'N', b'3', b' '),
    ),
    (
        bcnam(b' ', b'D', b'C'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b'D', b'C'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b' ', b'T'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b' ', b'T'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'O', b'4', b' '),
    ),
    (
        bcnam(b' ', b' ', b'T'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b' ', b'T'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b'D', b'T'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b'D', b'T'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'O', b'4', b' '),
    ),
    (
        bcnam(b' ', b'D', b'T'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b'D', b'T'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b' ', b'U'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b' ', b'U'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'O', b'4', b' '),
    ),
    (
        bcnam(b' ', b' ', b'U'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b' ', b'U'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
    (
        bcnam(b' ', b'D', b'U'),
        bcatm(b' ', b'C', b'2', b' '),
        bcatm(b' ', b'O', b'2', b' '),
    ),
    (
        bcnam(b' ', b'D', b'U'),
        bcatm(b' ', b'C', b'4', b' '),
        bcatm(b' ', b'O', b'4', b' '),
    ),
    (
        bcnam(b' ', b'D', b'U'),
        bcatm(b' ', b'C', b'5', b' '),
        bcatm(b' ', b'C', b'6', b' '),
    ),
    (
        bcnam(b' ', b'D', b'U'),
        bcatm(b' ', b'O', b'P', b'1'),
        bcatm(b' ', b'P', b' ', b' '),
    ),
];

fn standard_pdb_double_bond_like_rdkit(
    builder: &MoleculeBuilder,
    begin: AtomId,
    end: AtomId,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION StandardPDBDoubleBond atom wrapper
    // RDKit✔️✔️: static bool StandardPDBDoubleBond(RWMol *mol, Atom *beg, Atom *end) {
    // RDKit✔️✔️:   auto *bInfo = (AtomPDBResidueInfo *)beg->getMonomerInfo();
    // RDKit✔️✔️:   if (!bInfo || bInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto *eInfo = (AtomPDBResidueInfo *)end->getMonomerInfo();
    // RDKit✔️✔️:   if (!eInfo || eInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!SamePDBResidue(bInfo, eInfo)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bInfo->getIsHeteroAtom() || eInfo->getIsHeteroAtom()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const char *ptr = bInfo->getResidueName().c_str();
    // RDKit✔️✔️:   unsigned int rescode = BCNAM(ptr[0], ptr[1], ptr[2]);
    // RDKit✔️✔️:   ptr = bInfo->getName().c_str();
    // RDKit✔️✔️:   unsigned int atm1 = BCATM(ptr[0], ptr[1], ptr[2], ptr[3]);
    // RDKit✔️✔️:   ptr = eInfo->getName().c_str();
    // RDKit✔️✔️:   unsigned int atm2 = BCATM(ptr[0], ptr[1], ptr[2], ptr[3]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!StandardPDBDoubleBond(rescode, atm1, atm2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Check that neither end already has a double bond
    // RDKit✔️✔️:   ROMol::OBOND_ITER_PAIR bp;
    // RDKit✔️✔️:   for (bp = mol->getAtomBonds(beg); bp.first != bp.second; ++bp.first) {
    // RDKit✔️✔️:     if ((*mol)[*bp.first]->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (bp = mol->getAtomBonds(end); bp.first != bp.second; ++bp.first) {
    // RDKit✔️✔️:     if ((*mol)[*bp.first]->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StandardPDBDoubleBond atom wrapper
    let Some(begin_atom) = builder.atoms().get(begin.index()) else {
        return false;
    };
    let Some(end_atom) = builder.atoms().get(end.index()) else {
        return false;
    };
    let Some(begin_info) = begin_atom.pdb_residue_info() else {
        return false;
    };
    let Some(end_info) = end_atom.pdb_residue_info() else {
        return false;
    };
    if !same_pdb_residue_like_rdkit(begin_info, end_info) {
        return false;
    }
    if begin_info.is_hetero_atom() || end_info.is_hetero_atom() {
        return false;
    }
    let rescode = fixed_width_code(begin_info.residue_name().as_bytes(), 3);
    let atm1 = fixed_width_code(begin_info.atom_name().as_bytes(), 4);
    let atm2 = fixed_width_code(end_info.atom_name().as_bytes(), 4);
    if !standard_pdb_double_bond_table_like_rdkit(rescode, atm1, atm2) {
        return false;
    }
    for bond_id in builder.neighbor_bonds(begin) {
        if builder
            .bond(*bond_id)
            .is_some_and(|bond| bond.order() == BondOrder::Double)
        {
            return false;
        }
    }
    for bond_id in builder.neighbor_bonds(end) {
        if builder
            .bond(*bond_id)
            .is_some_and(|bond| bond.order() == BondOrder::Double)
        {
            return false;
        }
    }
    true
}

fn standard_pdb_residue_bond_orders_like_rdkit(
    builder: &mut MoleculeBuilder,
) -> Result<(), MoleculeBuildError> {
    // BEGIN RDKIT CPP FUNCTION StandardPDBResidueBondOrders
    // RDKit✔️✔️: void StandardPDBResidueBondOrders(RWMol *mol) {
    // RDKit✔️✔️:   RWMol::BondIterator bondIt;
    // RDKit✔️✔️:   for (bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    // RDKit✔️✔️:     Bond *bond = *bondIt;
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       Atom *beg = bond->getBeginAtom();
    // RDKit✔️✔️:       Atom *end = bond->getEndAtom();
    // RDKit✔️✔️:       if (StandardPDBDoubleBond(mol, beg, end)) {
    // RDKit✔️✔️:         bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StandardPDBResidueBondOrders
    for bond_index in 0..builder.bonds().len() {
        let bond_id = crate::BondId::new(bond_index);
        let Some(bond) = builder.bond(bond_id) else {
            continue;
        };
        if bond.order() != BondOrder::Single {
            continue;
        }
        if standard_pdb_double_bond_like_rdkit(builder, bond.begin(), bond.end()) {
            builder.set_bond_order(bond_id, BondOrder::Double)?;
        }
    }
    Ok(())
}

fn basic_pdb_cleanup_like_rdkit(builder: &mut MoleculeBuilder) {
    // BEGIN RDKIT CPP FUNCTION BasicPDBCleanup
    // RDKit✔️✔️: void BasicPDBCleanup(RWMol &mol) {
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    // RDKit✔️✔️:     Atom *atom = mol[*atBegin];
    // RDKit❗✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 &&
    // RDKit✔️✔️:         atom->getValence(Atom::ValenceType::EXPLICIT) == 4) {
    // RDKit✔️✔️:       atom->setFormalCharge(1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atBegin;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BasicPDBCleanup
    let atom_count = builder.atoms().len();
    for atom_index in 0..atom_count {
        let atom_id = AtomId::new(atom_index);
        let atom = &builder.atoms()[atom_index];
        if atom.atomic_number() != 7 || atom.formal_charge() != 0 {
            continue;
        }
        let explicit_valence = builder
            .neighbor_bonds(atom_id)
            .iter()
            .filter_map(|bond| builder.bond(*bond))
            .map(|bond| match bond.order() {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
                _ => 0,
            })
            .sum::<i32>();
        if explicit_valence == 4
            && let Some(atom) = builder.atom_mut(atom_id)
        {
            atom.set_formal_charge(1);
        }
    }
}

fn standard_pdb_chiral_atom_like_rdkit(residue_name: &str, atom_name: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION StandardPDBChiralAtom
    // RDKit✔️✔️: bool StandardPDBChiralAtom(const char *resnam, const char *atmnam) {
    // RDKit✔️✔️:   switch (BCNAM(resnam[0], resnam[1], resnam[2])) {
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'Y'):
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     case BCNAM('I', 'L', 'E'):
    // RDKit✔️✔️:     case BCNAM('T', 'H', 'R'):
    // RDKit✔️✔️:       // Alpha and beta carbons (" CA " and " CB ").
    // RDKit✔️✔️:       return atmnam[0] == ' ' && atmnam[1] == 'C' &&
    // RDKit✔️✔️:              (atmnam[2] == 'A' || atmnam[2] == 'B') && atmnam[3] == ' ';
    // RDKit✔️✔️:     case BCNAM('A', 'L', 'A'):
    // RDKit✔️✔️:     case BCNAM('A', 'R', 'G'):
    // RDKit✔️✔️:     case BCNAM('A', 'S', 'N'):
    // RDKit✔️✔️:     case BCNAM('A', 'S', 'P'):
    // RDKit✔️✔️:     case BCNAM('C', 'Y', 'S'):
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'N'):
    // RDKit✔️✔️:     case BCNAM('G', 'L', 'U'):
    // RDKit✔️✔️:     case BCNAM('H', 'I', 'S'):
    // RDKit✔️✔️:     case BCNAM('L', 'E', 'U'):
    // RDKit✔️✔️:     case BCNAM('L', 'Y', 'S'):
    // RDKit✔️✔️:     case BCNAM('M', 'E', 'T'):
    // RDKit✔️✔️:     case BCNAM('P', 'H', 'E'):
    // RDKit✔️✔️:     case BCNAM('P', 'R', 'O'):
    // RDKit✔️✔️:     case BCNAM('S', 'E', 'R'):
    // RDKit✔️✔️:     case BCNAM('T', 'R', 'P'):
    // RDKit✔️✔️:     case BCNAM('T', 'Y', 'R'):
    // RDKit✔️✔️:     case BCNAM('V', 'A', 'L'):
    // RDKit✔️✔️:       return atmnam[0] == ' ' && atmnam[1] == 'C' && atmnam[2] == 'A' &&
    // RDKit✔️✔️:              atmnam[3] == ' ';
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StandardPDBChiralAtom
    let rescode = fixed_width_code(residue_name.as_bytes(), 3);
    let atom = atom_name.as_bytes();
    match rescode {
        res if res == bcnam(b'G', b'L', b'Y') => false,
        res if res == bcnam(b'I', b'L', b'E') || res == bcnam(b'T', b'H', b'R') => {
            atom.get(0) == Some(&b' ')
                && atom.get(1) == Some(&b'C')
                && matches!(atom.get(2), Some(b'A' | b'B'))
                && atom.get(3) == Some(&b' ')
        }
        res if matches!(
            res,
            x if x == bcnam(b'A', b'L', b'A')
                || x == bcnam(b'A', b'R', b'G')
                || x == bcnam(b'A', b'S', b'N')
                || x == bcnam(b'A', b'S', b'P')
                || x == bcnam(b'C', b'Y', b'S')
                || x == bcnam(b'G', b'L', b'N')
                || x == bcnam(b'G', b'L', b'U')
                || x == bcnam(b'H', b'I', b'S')
                || x == bcnam(b'L', b'E', b'U')
                || x == bcnam(b'L', b'Y', b'S')
                || x == bcnam(b'M', b'E', b'T')
                || x == bcnam(b'P', b'H', b'E')
                || x == bcnam(b'P', b'R', b'O')
                || x == bcnam(b'S', b'E', b'R')
                || x == bcnam(b'T', b'R', b'P')
                || x == bcnam(b'T', b'Y', b'R')
                || x == bcnam(b'V', b'A', b'L')
        ) =>
        {
            atom.get(0) == Some(&b' ')
                && atom.get(1) == Some(&b'C')
                && atom.get(2) == Some(&b'A')
                && atom.get(3) == Some(&b' ')
        }
        _ => false,
    }
}

fn standard_pdb_residue_chirality_like_rdkit(molecule: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION StandardPDBResidueChirality
    // RDKit✔️✔️: void StandardPDBResidueChirality(RWMol *mol) {
    // RDKit✔️✔️:   for (ROMol::AtomIterator atomIt = mol->beginAtoms();
    // RDKit✔️✔️:        atomIt != mol->endAtoms(); ++atomIt) {
    // RDKit✔️✔️:     Atom *atom = *atomIt;
    // RDKit✔️✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:       auto *info = (AtomPDBResidueInfo *)atom->getMonomerInfo();
    // RDKit✔️✔️:       if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
    // RDKit✔️✔️:           !info->getIsHeteroAtom() &&
    // RDKit✔️✔️:           !StandardPDBChiralAtom(info->getResidueName().c_str(),
    // RDKit✔️✔️:                                  info->getName().c_str())) {
    // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:         if (atom->hasProp(common_properties::_CIPCode)) {
    // RDKit✔️✔️:           atom->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StandardPDBResidueChirality
    for atom in &mut molecule.topology_block_mut().atoms {
        if atom.chiral_tag() == ChiralTag::Unspecified {
            continue;
        }
        let should_clear = atom.pdb_residue_info().is_some_and(|info| {
            !info.is_hetero_atom()
                && !standard_pdb_chiral_atom_like_rdkit(info.residue_name(), info.atom_name())
        });
        if should_clear {
            atom.set_chiral_tag(ChiralTag::Unspecified);
            atom.clear_prop("_CIPCode");
        }
    }
}

fn atom_name_string(name: crate::AtomName) -> String {
    String::from_utf8_lossy(&name.0).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_proximity_profile() -> RdkitPdbMolProfile {
        RdkitPdbMolProfile {
            sanitize: false,
            remove_hs: false,
            flavor: 0,
            proximity_bonding: false,
        }
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_retains_atom_and_hetatm_records() {
        let pdb = "\
ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N
HETATM    2  O   HOH B   2      12.000  14.000   8.000  0.50 10.00           O
";
        let mol = Molecule::from_pdb_block_with_options(pdb, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.atoms()[0].atomic_number(), 7);
        assert_eq!(mol.atoms()[1].atomic_number(), 8);
        let first = mol.atoms()[0].pdb_residue_info().unwrap();
        let second = mol.atoms()[1].pdb_residue_info().unwrap();
        assert!(!first.is_hetero_atom());
        assert!(second.is_hetero_atom());
        assert_eq!(second.residue_name(), "HOH");
        assert_eq!(second.occupancy(), 0.5);
        assert_eq!(second.temp_factor(), 10.0);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_filters_altloc_like_rdkit_default_flavor() {
        let pdb = "\
ATOM      1  CA AALA A   1      11.104  13.207   9.900  1.00 20.00           C
ATOM      2  CA BALA A   1      12.104  13.207   9.900  1.00 20.00           C
ATOM      3  CA 1ALA A   1      13.104  13.207   9.900  1.00 20.00           C
";
        let mol = Molecule::from_pdb_block_with_options(pdb, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(
            mol.atoms()
                .iter()
                .map(|atom| atom.pdb_residue_info().unwrap().serial_number())
                .collect::<Vec<_>>(),
            vec![1, 3]
        );
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_flavor_one_keeps_altlocs() {
        let pdb = "\
ATOM      1  CA AALA A   1      11.104  13.207   9.900  1.00 20.00           C
ATOM      2  CA BALA A   1      12.104  13.207   9.900  1.00 20.00           C
";
        let mol = Molecule::from_pdb_block_with_options(
            pdb,
            RdkitPdbMolProfile {
                flavor: 1,
                ..no_proximity_profile()
            },
        )
        .unwrap();

        assert_eq!(mol.num_atoms(), 2);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_applies_conect_duplicate_bond_orders() {
        let pdb = "\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C
HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00 10.00           O
CONECT    1    2    2
";
        let mol = Molecule::from_pdb_block_with_options(pdb, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bonds()[0].order(), BondOrder::Double);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_adds_proximity_bonds_like_rdkit() {
        let pdb = "\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C
HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00 10.00           O
";
        let mol = Molecule::from_pdb_block_with_options(
            pdb,
            RdkitPdbMolProfile {
                sanitize: false,
                remove_hs: false,
                flavor: 0,
                proximity_bonding: true,
            },
        )
        .unwrap();

        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bonds()[0].order(), BondOrder::Single);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_default_proximity_ignores_h_h_contacts() {
        let pdb = "\
HETATM    1  H1  LIG A   1       0.000   0.000   0.000  1.00 10.00           H
HETATM    2  H2  LIG A   1       0.700   0.000   0.000  1.00 10.00           H
";
        let mol = Molecule::from_pdb_block_with_options(
            pdb,
            RdkitPdbMolProfile {
                sanitize: false,
                remove_hs: false,
                flavor: 0,
                proximity_bonding: true,
            },
        )
        .unwrap();

        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_blacklisted_explicit_conect_becomes_zero_bond() {
        let pdb = "\
HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00 10.00          FE
HETATM    2  O   HOH B   2       2.000   0.000   0.000  1.00 10.00           O
CONECT    1    2
";
        let mol = Molecule::from_pdb_block_with_options(pdb, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bonds()[0].order(), BondOrder::Zero);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_accepts_full_pdb_element_table() {
        let pdb = "\
HETATM    1 HG    HG A   1      -2.213  10.563  24.265  1.00 32.73          HG
HETATM    2 CD    CD A   2      -3.467  18.396  77.649  0.50 39.48          CD
";
        let mol = Molecule::from_pdb_block_with_options(pdb, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.atoms()[0].atomic_number(), 80);
        assert_eq!(mol.atoms()[1].atomic_number(), 48);
    }

    #[test]
    fn rdkit_pdb_molecule_conversion_flavor_eight_applies_standard_residue_double_bonds() {
        let pdb = "\
ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  O   ALA A   1       1.200   0.000   0.000  1.00 10.00           O
CONECT    1    2
";
        let mol = Molecule::from_pdb_block_with_options(
            pdb,
            RdkitPdbMolProfile {
                flavor: 8,
                ..no_proximity_profile()
            },
        )
        .unwrap();

        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bonds()[0].order(), BondOrder::Double);
    }

    #[test]
    fn mmcif_molecule_conversion_uses_bio_structure_then_rdkit_pdb_profile() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
ATOM 1 C CA . ALA A 1 7 ? 11.104 13.207 9.900 0.50 10.00 7 ALA A CA
HETATM 2 O O . HOH B 2 2 ? 12.000 14.000 8.000 1.00 20.00 2 HOH B O
"#;

        let mol = Molecule::from_mmcif_block_with_options(cif, no_proximity_profile()).unwrap();

        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 0);
        assert_eq!(mol.atoms()[0].atomic_number(), 6);
        assert_eq!(mol.atoms()[1].atomic_number(), 8);
        let first = mol.atoms()[0].pdb_residue_info().unwrap();
        let second = mol.atoms()[1].pdb_residue_info().unwrap();
        assert_eq!(first.atom_name(), "CA  ");
        assert_eq!(first.serial_number(), 1);
        assert_eq!(first.residue_name(), "ALA");
        assert_eq!(first.residue_number(), 7);
        assert_eq!(first.chain_id(), "A");
        assert_eq!(first.insertion_code(), " ");
        assert_eq!(first.occupancy(), 0.5);
        assert_eq!(first.temp_factor(), 10.0);
        assert!(!first.is_hetero_atom());
        assert_eq!(second.residue_name(), "HOH");
        assert_eq!(second.chain_id(), "B");
        assert!(second.is_hetero_atom());
        let coords = mol.conformers_3d()[0].coordinates();
        for (actual, expected) in coords
            .iter()
            .zip([[11.104, 13.207, 9.9], [12.0, 14.0, 8.0]])
        {
            for (actual_component, expected_component) in actual.iter().zip(expected) {
                assert!((actual_component - expected_component).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn mmcif_molecule_conversion_applies_same_proximity_profile() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 . LIG A 1 0.000 0.000 0.000
HETATM 2 O O1 . LIG A 1 1.200 0.000 0.000
"#;

        let mol = Molecule::from_mmcif_block_with_options(
            cif,
            RdkitPdbMolProfile {
                sanitize: false,
                remove_hs: false,
                flavor: 0,
                proximity_bonding: true,
            },
        )
        .unwrap();

        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bonds()[0].order(), BondOrder::Single);
    }
}
