// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// RDKit source: third_party/rdkit/Code/GraphMol/FileParsers/PDBWriter.cpp
// RDKit source: third_party/rdkit/Code/GraphMol/FileParsers/PDBParser.cpp
//
// Port of PDBWriter and MolToPDBBlock from RDKit to COSMolKit.
//
// PDBWriter support multiple "flavors" of PDB output
// flavor & 1 : Write MODEL/ENDMDL lines around each record
// flavor & 2 : Don't write any CONECT records
// flavor & 4 : Write CONECT records in both directions
// flavor & 8 : Don't use multiple CONECTs to encode bond order
// flavor & 16 : Write MASTER record
// flavor & 32 : Write TER record

use std::collections::BTreeMap;

use crate::{Atom, AtomId, BondOrder, Conformer3D, Molecule};

// RDKit source: PDBWriter.cpp lines 42-43
// RDKit✔️✔️: std::string GetDefaultAtomNumber(const Atom *atom,
// RDKit✔️✔️:                                  std::map<unsigned int, unsigned int> &elem);

// RDKit source: PDBWriter.cpp lines 141-165
// RDKit✔️✔️: std::string GetDefaultAtomNumber(const Atom *atom,
// RDKit✔️✔️:                                  std::map<unsigned int, unsigned int> &elem) {
// RDKit✔️✔️:   std::string ret = "  ";
// RDKit✔️✔️:   unsigned int atno = atom->getAtomicNum();
// RDKit✔️✔️:   if (elem.find(atno) == elem.end()) {
// RDKit✔️✔️:     elem[atno] = 1;
// RDKit✔️✔️:     ret[0] = '1';
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     unsigned int tmp = elem[atno] + 1;
// RDKit✔️✔️:     elem[atno] = tmp;
// RDKit✔️✔️:     if (tmp < 10) {
// RDKit✔️✔️:       ret[0] = tmp + '0';
// RDKit✔️✔️:     } else if (tmp < 100) {
// RDKit✔️✔️:       ret[0] = (tmp / 10) + '0';
// RDKit✔️✔️:       ret[1] = (tmp % 10) + '0';
// RDKit✔️✔️:     } else if (tmp < 360) {
// RDKit✔️✔️:       ret[0] = ((tmp - 100) / 10) + 'A';
// RDKit✔️✔️:       ret[1] = ((tmp - 100) % 10) + '0';
// RDKit✔️✔️:     } else if (tmp < 1036) {
// RDKit✔️✔️:       ret[0] = ((tmp - 360) / 26) + 'A';
// RDKit✔️✔️:       ret[1] = ((tmp - 360) % 26) + 'A';
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return ret;
// RDKit✔️✔️: }
//
// Track element counts and format a 2-character atom number string.
fn get_default_atom_number(atomic_number: u8, elem: &mut BTreeMap<u8, u32>) -> String {
    let mut ret = [' ', ' '];
    let entry = elem.entry(atomic_number).or_insert(0);
    *entry += 1;
    let tmp = *entry;
    if tmp < 10 {
        ret[0] = char::from_u32(u32::from(b'0') + tmp).unwrap();
    } else if tmp < 100 {
        ret[0] = char::from_u32(u32::from(b'0') + tmp / 10).unwrap();
        ret[1] = char::from_u32(u32::from(b'0') + tmp % 10).unwrap();
    } else if tmp < 360 {
        ret[0] = char::from_u32(u32::from(b'A') + (tmp - 100) / 10).unwrap();
        ret[1] = char::from_u32(u32::from(b'0') + (tmp - 100) % 10).unwrap();
    } else if tmp < 1036 {
        ret[0] = char::from_u32(u32::from(b'A') + (tmp - 360) / 26).unwrap();
        ret[1] = char::from_u32(u32::from(b'A') + (tmp - 360) % 26).unwrap();
    }
    ret.iter().collect()
}

/// Map atomic number to element symbol.
/// Covers elements 0-118 (H through Og). Unknown numbers return "?".
fn element_symbol(atomic_num: u8) -> &'static str {
    // COSMolKit's Element currently permits u8 values outside RDKit's table.
    // Preserve this writer's existing fallback at the format boundary.
    crate::rdkit_element_symbol(atomic_num).unwrap_or("?")
}
// RDKit source: PDBWriter.cpp lines 45-139
// RDKit✔️✔️: std::string GetPDBAtomLine(const Atom *atom, const Conformer *conf,
// RDKit✔️✔️:                            std::map<unsigned int, unsigned int> &elem) {
// RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
// RDKit✔️✔️:   std::stringstream ss;
// RDKit✔️✔️:
// RDKit✔️✔️:   std::string symb = atom->getSymbol();
// RDKit✔️✔️:   char at1, at2;
// RDKit✔️✔️:   switch (symb.length()) {
// RDKit✔️✔️:     case 0:
// RDKit✔️✔️:       at1 = ' ';
// RDKit✔️✔️:       at2 = 'X';
// RDKit✔️✔️:       break;
// RDKit✔️✔️:     case 1:
// RDKit✔️✔️:       at1 = ' ';
// RDKit✔️✔️:       at2 = symb[0];
// RDKit✔️✔️:       break;
// RDKit✔️✔️:     default:
// RDKit✔️✔️:       at1 = symb[0];
// RDKit✔️✔️:       at2 = symb[1];
// RDKit✔️✔️:       if (at2 >= 'a' && at2 <= 'z') {
// RDKit✔️✔️:         at2 -= 32;  // toupper
// RDKit✔️✔️:       }
// RDKit✔️✔️:       break;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   auto *info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
// RDKit✔️✔️:   if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
// RDKit✔️✔️:     ss << (info->getIsHeteroAtom() ? "HETATM" : "ATOM  ");
// RDKit✔️✔️:     ss << std::setw(5) << atom->getIdx() + 1;
// RDKit✔️✔️:     ss << ' ';
// RDKit✔️✔️:     const std::string &name = info->getName();
// RDKit✔️✔️:     if (name.empty()) {
// RDKit✔️✔️:       std::string atnum = GetDefaultAtomNumber(atom, elem);
// RDKit✔️✔️:       ss << at1 << at2 << atnum;
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       ss << std::setw(4) << name.substr(0, 4);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     const char *ptr = info->getAltLoc().c_str();
// RDKit✔️✔️:     if (*ptr == '\0') {
// RDKit✔️✔️:       ptr = " ";
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ss << *ptr;
// RDKit✔️✔️:     ss << std::setw(3) << info->getResidueName().substr(0, 3);
// RDKit✔️✔️:     ss << ' ';
// RDKit✔️✔️:     ptr = info->getChainId().c_str();
// RDKit✔️✔️:     if (*ptr == '\0') {
// RDKit✔️✔️:       ptr = " ";
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ss << *ptr;
// RDKit✔️✔️:     ss << std::setw(4) << info->getResidueNumber();
// RDKit✔️✔️:     ptr = info->getInsertionCode().c_str();
// RDKit✔️✔️:     if (*ptr == '\0') {
// RDKit✔️✔️:       ptr = " ";
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ss << *ptr;
// RDKit✔️✔️:     ss << "   ";
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     info = (AtomPDBResidueInfo *)nullptr;
// RDKit✔️✔️:     std::string atnum = GetDefaultAtomNumber(atom, elem);
// RDKit✔️✔️:     ss << "HETATM";
// RDKit✔️✔️:     ss << std::setw(5) << atom->getIdx() + 1;
// RDKit✔️✔️:     ss << ' ';
// RDKit✔️✔️:     ss << at1 << at2 << atnum;
// RDKit✔️✔️:     ss << " UNL     1    ";
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   if (conf) {
// RDKit✔️✔️:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
// RDKit✔️✔️:     ss << boost::format("%8.3f%8.3f%8.3f") % pos.x % pos.y % pos.z;
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     ss << "   0.000   0.000   0.000";
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   if (info) {
// RDKit✔️✔️:     ss << boost::format("%6.2f%6.2f") % info->getOccupancy() %
// RDKit✔️✔️:               info->getTempFactor();
// RDKit✔️✔️:     ss << "          ";
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     ss << "  1.00  0.00          ";
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   ss << at1;
// RDKit✔️✔️:   ss << at2;
// RDKit✔️✔️:   int charge = atom->getFormalCharge();
// RDKit✔️✔️:   if (charge > 0 && charge < 10) {
// RDKit✔️✔️:     ss << (char)('0' + charge);
// RDKit✔️✔️:     ss << '+';
// RDKit✔️✔️:   } else if (charge < 0 && charge > -10) {
// RDKit✔️✔️:     ss << (char)('0' - charge);
// RDKit✔️✔️:     ss << '-';
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     ss << "  ";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return ss.str();
// RDKit✔️✔️: }
//
// Build a single PDB ATOM/HETATM line for the given atom.
fn get_pdb_atom_line(atom: &Atom, conf: Option<&Conformer3D>, elem: &mut BTreeMap<u8, u32>) -> String {
    let symb = element_symbol(atom.atomic_number());
    let (at1, at2): (char, char) = match symb.len() {
        0 => (' ', 'X'),
        1 => (' ', symb.chars().next().unwrap()),
        _ => {
            let chars: Vec<char> = symb.chars().collect();
            let first = chars[0];
            let mut second = chars[1];
            if second.is_ascii_lowercase() {
                second = second.to_ascii_uppercase();
            }
            (first, second)
        }
    };

    let info = atom.pdb_residue_info();
    let mut line = String::with_capacity(80);

    if let Some(pdb_info) = info {
        // RDKit✔️✔️: ss << (info->getIsHeteroAtom() ? "HETATM" : "ATOM  ");
        if pdb_info.is_hetero_atom() {
            line.push_str("HETATM");
        } else {
            line.push_str("ATOM  ");
        }
        // RDKit✔️✔️: ss << std::setw(5) << atom->getIdx() + 1;
        line.push_str(&format!("{:>5}", atom.id().index() + 1));
        line.push(' ');
        // RDKit✔️✔️: const std::string &name = info->getName();
        let name = pdb_info.atom_name();
        if name.is_empty() {
            // RDKit✔️✔️: std::string atnum = GetDefaultAtomNumber(atom, elem);
            // RDKit✔️✔️: ss << at1 << at2 << atnum;
            let atnum = get_default_atom_number(atom.atomic_number(), elem);
            line.push(at1);
            line.push(at2);
            line.push_str(&atnum);
        } else {
            // RDKit✔️✔️: ss << std::setw(4) << name.substr(0, 4);
            let truncated: String = name.chars().take(4).collect();
            line.push_str(&format!("{:>4}", truncated));
        }
        // RDKit✔️✔️: const char *ptr = info->getAltLoc().c_str();
        // RDKit✔️✔️: if (*ptr == '\0') { ptr = " "; }
        // RDKit✔️✔️: ss << *ptr;
        let alt = pdb_info.alt_loc();
        line.push(alt.chars().next().unwrap_or(' '));
        // RDKit✔️✔️: ss << std::setw(3) << info->getResidueName().substr(0, 3);
        let resname: String = pdb_info.residue_name().chars().take(3).collect();
        line.push_str(&format!("{:>3}", resname));
        line.push(' ');
        // RDKit✔️✔️: ptr = info->getChainId().c_str();
        // RDKit✔️✔️: if (*ptr == '\0') { ptr = " "; }
        // RDKit✔️✔️: ss << *ptr;
        let chain = pdb_info.chain_id();
        line.push(chain.chars().next().unwrap_or(' '));
        // RDKit✔️✔️: ss << std::setw(4) << info->getResidueNumber();
        line.push_str(&format!("{:>4}", pdb_info.residue_number()));
        // RDKit✔️✔️: ptr = info->getInsertionCode().c_str();
        // RDKit✔️✔️: if (*ptr == '\0') { ptr = " "; }
        // RDKit✔️✔️: ss << *ptr;
        let icode = pdb_info.insertion_code();
        line.push(icode.chars().next().unwrap_or(' '));
        line.push_str("   ");
    } else {
        // RDKit✔️✔️: std::string atnum = GetDefaultAtomNumber(atom, elem);
        let atnum = get_default_atom_number(atom.atomic_number(), elem);
        // RDKit✔️✔️: ss << "HETATM";
        line.push_str("HETATM");
        // RDKit✔️✔️: ss << std::setw(5) << atom->getIdx() + 1;
        line.push_str(&format!("{:>5}", atom.id().index() + 1));
        line.push(' ');
        // RDKit✔️✔️: ss << at1 << at2 << atnum;
        line.push(at1);
        line.push(at2);
        line.push_str(&atnum);
        // RDKit✔️✔️: ss << " UNL     1    ";
        line.push_str(" UNL     1    ");
    }

    // RDKit✔️✔️: if (conf) {
    // RDKit✔️✔️:   const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit✔️✔️:   ss << boost::format("%8.3f%8.3f%8.3f") % pos.x % pos.y % pos.z;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   ss << "   0.000   0.000   0.000";
    // RDKit✔️✔️: }
    if let Some(conf3d) = conf {
        let coords = conf3d.coordinates();
        if atom.id().index() < coords.len() {
            let pos = coords[atom.id().index()];
            line.push_str(&format!("{:+8.3}{:+8.3}{:+8.3}", pos[0], pos[1], pos[2]));
        } else {
            line.push_str("   0.000   0.000   0.000");
        }
    } else {
        line.push_str("   0.000   0.000   0.000");
    }

    // RDKit✔️✔️: if (info) {
    // RDKit✔️✔️:   ss << boost::format("%6.2f%6.2f") % info->getOccupancy() %
    // RDKit✔️✔️:             info->getTempFactor();
    // RDKit✔️✔️:   ss << "          ";
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   ss << "  1.00  0.00          ";
    // RDKit✔️✔️: }
    if let Some(pdb_info) = info {
        line.push_str(&format!("{:>6.2}{:>6.2}", pdb_info.occupancy(), pdb_info.temp_factor()));
        line.push_str("          ");
    } else {
        line.push_str("  1.00  0.00          ");
    }

    // RDKit✔️✔️: ss << at1;
    // RDKit✔️✔️: ss << at2;
    line.push(at1);
    line.push(at2);
    // RDKit✔️✔️: int charge = atom->getFormalCharge();
    // RDKit✔️✔️: if (charge > 0 && charge < 10) {
    // RDKit✔️✔️:   ss << (char)('0' + charge);
    // RDKit✔️✔️:   ss << '+';
    // RDKit✔️✔️: } else if (charge < 0 && charge > -10) {
    // RDKit✔️✔️:   ss << (char)('0' - charge);
    // RDKit✔️✔️:   ss << '-';
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   ss << "  ";
    // RDKit✔️✔️: }
    let charge = atom.formal_charge();
    if charge > 0 && charge < 10 {
        line.push(char::from_u32(u32::from(b'0') + charge as u32).unwrap());
        line.push('+');
    } else if charge < 0 && charge > -10 {
        line.push(char::from_u32(u32::from(b'0') + (-charge) as u32).unwrap());
        line.push('-');
    } else {
        line.push(' ');
        line.push(' ');
    }

    line
}

// RDKit source: PDBWriter.cpp lines 167-226
// RDKit✔️✔️: std::string GetPDBBondLines(const Atom *atom, bool all, bool both, bool mult,
// RDKit✔️✔️:                             unsigned int &conect_count) {
// RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
// RDKit✔️✔️:   unsigned int src = atom->getIdx() + 1;
// RDKit✔️✔️:   std::vector<unsigned int> v;
// RDKit✔️✔️:
// RDKit✔️✔️:   ROMol *mol = &atom->getOwningMol();
// RDKit✔️✔️:   for (ROMol::OBOND_ITER_PAIR bondIt = mol->getAtomBonds(atom);
// RDKit✔️✔️:        bondIt.first != bondIt.second; ++bondIt.first) {
// RDKit✔️✔️:     Bond *bptr = (*mol)[*bondIt.first];
// RDKit✔️✔️:     Atom *nptr = bptr->getOtherAtom(atom);
// RDKit✔️✔️:     unsigned int dst = nptr->getIdx() + 1;
// RDKit✔️✔️:     if (dst < src && !both) {
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     Bond::BondType btype = Bond::SINGLE;
// RDKit✔️✔️:     if (mult) {
// RDKit✔️✔️:       btype = bptr->getBondType();
// RDKit✔️✔️:     }
// RDKit✔️✔️:     switch (btype) {
// RDKit✔️✔️:       default:
// RDKit✔️✔️:       case Bond::SINGLE:
// RDKit✔️✔️:       case Bond::AROMATIC:
// RDKit✔️✔️:         if (all) {
// RDKit✔️✔️:           v.push_back(dst);
// RDKit✔️✔️:         }
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       case Bond::QUADRUPLE:
// RDKit✔️✔️:         v.push_back(dst);
// RDKit✔️✔️:         /* FALLTHRU */
// RDKit✔️✔️:       case Bond::TRIPLE:
// RDKit✔️✔️:         v.push_back(dst);
// RDKit✔️✔️:         /* FALLTHRU */
// RDKit✔️✔️:       case Bond::DOUBLE:
// RDKit✔️✔️:         v.push_back(dst);
// RDKit✔️✔️:         v.push_back(dst);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   auto count = rdcast<unsigned int>(v.size());
// RDKit✔️✔️:   if (count == 0) {
// RDKit✔️✔️:     return "";
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   std::sort(v.begin(), v.end());
// RDKit✔️✔️:   std::stringstream ss;
// RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
// RDKit✔️✔️:     if ((i & 3) == 0) {
// RDKit✔️✔️:       if (i != 0) {
// RDKit✔️✔️:         ss << '\n';
// RDKit✔️✔️:       }
// RDKit✔️✔️:       ss << "CONECT";
// RDKit✔️✔️:       ss << std::setw(5) << src;
// RDKit✔️✔️:       conect_count++;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ss << std::setw(5) << v[i];
// RDKit✔️✔️:   }
// RDKit✔️✔️:   ss << '\n';
// RDKit✔️✔️:   return ss.str();
// RDKit✔️✔️: }
//
// Build CONECT lines for a single atom's bonds.
fn get_pdb_bond_lines(
    atom_idx: AtomId,
    bonds: &[crate::Bond],
    all: bool,
    both: bool,
    mult: bool,
    conect_count: &mut u32,
) -> String {
    let src = atom_idx.index() + 1;
    let mut v: Vec<usize> = Vec::new();

    // Scan bonds for those attached to this atom.
    for bond in bonds {
        let begin = bond.begin();
        let end = bond.end();
        let other_idx = if begin == atom_idx {
            end
        } else if end == atom_idx {
            begin
        } else {
            continue;
        };
        let dst = other_idx.index() + 1;

        // RDKit✔️✔️: if (dst < src && !both) { continue; }
        if dst < src && !both {
            continue;
        }

        // RDKit✔️✔️: Bond::BondType btype = Bond::SINGLE;
        // RDKit✔️✔️: if (mult) {
        // RDKit✔️✔️:   btype = bptr->getBondType();
        // RDKit✔️✔️: }
        let btype = if mult { bond.order() } else { BondOrder::Single };

        // RDKit✔️✔️: switch (btype) {
        // RDKit✔️✔️:   default:
        // RDKit✔️✔️:   case Bond::SINGLE:
        // RDKit✔️✔️:   case Bond::AROMATIC:
        // RDKit✔️✔️:     if (all) { v.push_back(dst); }
        // RDKit✔️✔️:     break;
        // RDKit✔️✔️:   case Bond::QUADRUPLE:  v.push_back(dst); /* FALLTHRU */
        // RDKit✔️✔️:   case Bond::TRIPLE:      v.push_back(dst); /* FALLTHRU */
        // RDKit✔️✔️:   case Bond::DOUBLE:      v.push_back(dst); v.push_back(dst);
        // RDKit✔️✔️: }
        match btype {
            BondOrder::Single | BondOrder::Aromatic => {
                if all {
                    v.push(dst);
                }
            }
            BondOrder::Quadruple => {
                v.push(dst);
                // fall through to Triple
                v.push(dst);
                // fall through to Double
                v.push(dst);
                v.push(dst);
            }
            BondOrder::Triple => {
                v.push(dst);
                // fall through to Double
                v.push(dst);
                v.push(dst);
            }
            BondOrder::Double => {
                v.push(dst);
                v.push(dst);
            }
            // default: no bonds emitted
            _ => {}
        }
    }

    if v.is_empty() {
        return String::new();
    }

    // RDKit✔️✔️: std::sort(v.begin(), v.end());
    v.sort_unstable();

    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     if ((i & 3) == 0) {
    // RDKit✔️✔️:       if (i != 0) { ss << '\n'; }
    // RDKit✔️✔️:       ss << "CONECT" << std::setw(5) << src;
    // RDKit✔️✔️:       conect_count++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ss << std::setw(5) << v[i];
    // RDKit✔️✔️:   }
    let mut result = String::new();
    for (i, dst) in v.iter().enumerate() {
        if (i & 3) == 0 {
            if i != 0 {
                result.push('\n');
            }
            result.push_str(&format!("CONECT{:>5}", src));
            *conect_count += 1;
        }
        result.push_str(&format!("{:>5}", dst));
    }
    result.push('\n');

    result
}

// RDKit source: PDBWriter.cpp lines 228-265
// RDKit✔️✔️: std::string MolToPDBBody(const ROMol &mol, const Conformer *conf,
// RDKit✔️✔️:                          unsigned int flavor, unsigned int &atm_count,
// RDKit✔️✔️:                          unsigned int &ter_count, unsigned int &conect_count) {
// RDKit✔️✔️:   std::string res;
// RDKit✔️✔️:   std::string last;
// RDKit✔️✔️:   std::map<unsigned int, unsigned int> elem;
// RDKit✔️✔️:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
// RDKit✔️✔️:        atomIt != mol.endAtoms(); ++atomIt) {
// RDKit✔️✔️:     last = GetPDBAtomLine(*atomIt, conf, elem);
// RDKit✔️✔️:     res += last;
// RDKit✔️✔️:     res += '\n';
// RDKit✔️✔️:     atm_count++;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   if (ter_count == 0 && atm_count && (flavor & 32)) {
// RDKit✔️✔️:     std::stringstream ss;
// RDKit✔️✔️:     ss << "TER   ";
// RDKit✔️✔️:     ss << std::setw(5) << atm_count + 1;
// RDKit✔️✔️:     if (last.length() >= 27) {
// RDKit✔️✔️:       ss << "      ";
// RDKit✔️✔️:       ss << last.substr(17, 10);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ss << '\n';
// RDKit✔️✔️:     res += ss.str();
// RDKit✔️✔️:     ter_count = 1;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   bool all = (flavor & 2) == 0;
// RDKit✔️✔️:   bool both = (flavor & 4) != 0;
// RDKit✔️✔️:   bool mult = (flavor & 8) == 0;
// RDKit✔️✔️:   if (all || mult) {
// RDKit✔️✔️:     for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
// RDKit✔️✔️:          atomIt != mol.endAtoms(); ++atomIt) {
// RDKit✔️✔️:       res += GetPDBBondLines(*atomIt, all, both, mult, conect_count);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
//
// Build the body of a PDB block: ATOM/HETATM + optional TER + CONECT records.
fn mol_to_pdb_body(
    mol: &Molecule,
    conf: Option<&Conformer3D>,
    flavor: u32,
    atm_count: &mut u32,
    ter_count: &mut u32,
    conect_count: &mut u32,
) -> String {
    let mut res = String::new();
    let mut last = String::new();
    let mut elem: BTreeMap<u8, u32> = BTreeMap::new();

    // RDKit✔️✔️: for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit✔️✔️:      atomIt != mol.endAtoms(); ++atomIt) {
    for atom in mol.atoms() {
        last = get_pdb_atom_line(atom, conf, &mut elem);
        res.push_str(&last);
        res.push('\n');
        *atm_count += 1;
    }

    // RDKit✔️✔️: if (ter_count == 0 && atm_count && (flavor & 32)) {
    if *ter_count == 0 && *atm_count > 0 && (flavor & 32) != 0 {
        let mut ter_line = String::with_capacity(80);
        // RDKit✔️✔️: ss << "TER   ";
        ter_line.push_str("TER   ");
        // RDKit✔️✔️: ss << std::setw(5) << atm_count + 1;
        ter_line.push_str(&format!("{:>5}", *atm_count + 1));
        // RDKit✔️✔️: if (last.length() >= 27) {
        // RDKit✔️✔️:   ss << "      ";
        // RDKit✔️✔️:   ss << last.substr(17, 10);
        // RDKit✔️✔️: }
        if last.len() >= 27 {
            // PDB columns 18-27 (0-indexed: 17..27) -> chain, resSeq, iCode, "   "
            // This is 10 chars starting at offset 17 in the 80-char PDB line.
            // In COSMolKit the atom line may be shorter; we take a substring.
            let end = (17 + 10).min(last.len());
            ter_line.push_str("      ");
            ter_line.push_str(&last[17..end]);
        }
        ter_line.push('\n');
        res.push_str(&ter_line);
        *ter_count = 1;
    }

    // RDKit✔️✔️: bool all = (flavor & 2) == 0;
    // RDKit✔️✔️: bool both = (flavor & 4) != 0;
    // RDKit✔️✔️: bool mult = (flavor & 8) == 0;
    let all = (flavor & 2) == 0;
    let both = (flavor & 4) != 0;
    let mult = (flavor & 8) == 0;

    // RDKit✔️✔️: if (all || mult) {
    if all || mult {
        for atom in mol.atoms() {
            res.push_str(&get_pdb_bond_lines(
                atom.id(),
                mol.bonds(),
                all,
                both,
                mult,
                conect_count,
            ));
        }
    }

    res
}

// RDKit source: PDBWriter.cpp lines 267-322
// RDKit✔️✔️: std::string MolToPDBBlock(const ROMol &imol, int confId, unsigned int flavor) {
// RDKit✔️✔️:   RWMol rwmol(imol);
// RDKit✔️✔️:   MolOps::Kekulize(rwmol);
// RDKit✔️✔️:   Utils::LocaleSwitcher ls;
// RDKit✔️✔️:
// RDKit✔️✔️:   std::string res;
// RDKit✔️✔️:   std::string name;
// RDKit✔️✔️:   if (rwmol.getPropIfPresent(common_properties::_Name, name)) {
// RDKit✔️✔️:     if (!name.empty()) {
// RDKit✔️✔️:       res += "COMPND    ";
// RDKit✔️✔️:       res += name;
// RDKit✔️✔️:       res += '\n';
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   unsigned int atm_count = 0;
// RDKit✔️✔️:   unsigned int ter_count = 0;
// RDKit✔️✔️:   unsigned int conect_count = 0;
// RDKit✔️✔️:
// RDKit✔️✔️:   const Conformer *conf;
// RDKit✔️✔️:   if (confId < 0 && rwmol.getNumConformers() > 1) {
// RDKit✔️✔️:     int count = rwmol.getNumConformers();
// RDKit✔️✔️:     for (confId = 0; confId < count; confId++) {
// RDKit✔️✔️:       conf = &(rwmol.getConformer(confId));
// RDKit✔️✔️:       std::stringstream ss;
// RDKit✔️✔️:       ss << "MODEL     ";
// RDKit✔️✔️:       ss << std::setw(4) << (confId + 1);
// RDKit✔️✔️:       ss << "\n";
// RDKit✔️✔️:       res += ss.str();
// RDKit✔️✔️:       res +=
// RDKit✔️✔️:           MolToPDBBody(rwmol, conf, flavor, atm_count, ter_count, conect_count);
// RDKit✔️✔️:       res += "ENDMDL\n";
// RDKit✔️✔️:     }
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     if (confId < 0 && rwmol.getNumConformers() == 0) {
// RDKit✔️✔️:       conf = nullptr;
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       conf = &(rwmol.getConformer(confId));
// RDKit✔️✔️:     }
// RDKit✔️✔️:     res +=
// RDKit✔️✔️:         MolToPDBBody(rwmol, conf, flavor, atm_count, ter_count, conect_count);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   if (flavor & 16) {
// RDKit✔️✔️:     std::stringstream ss;
// RDKit✔️✔️:     ss << "MASTER        0    0    0    0    0    0    0    0";
// RDKit✔️✔️:     ss << std::setw(5) << atm_count;
// RDKit✔️✔️:     ss << std::setw(5) << ter_count;
// RDKit✔️✔️:     ss << std::setw(5) << conect_count;
// RDKit✔️✔️:     ss << "    0\n";
// RDKit✔️✔️:     res += ss.str();
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   res += "END\n";
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
//
// Write a molecule to a PDB-format string block.
//
// `conf_id` selects the conformer index (0-based). Use -1 to write all
// conformers (with MODEL/ENDMDL wrappers) when multiple exist, or to use
// no coordinates when no conformers exist.
//
// `flavor` bitmask (see top-of-file comments):
//   1 = MODEL/ENDMDL per record, 2 = no CONECT, 4 = bidirectional CONECT,
//   8 = single CONECT per bond, 16 = MASTER record, 32 = TER record.
pub fn mol_to_pdb_block(mol: &Molecule, conf_id: i32, flavor: u32) -> String {
    // RDKit✔️✔️: std::string name;
    // RDKit✔️✔️: if (rwmol.getPropIfPresent(common_properties::_Name, name)) {
    // RDKit✔️✔️:   if (!name.empty()) {
    // RDKit✔️✔️:     res += "COMPND    ";
    // RDKit✔️✔️:     res += name;
    // RDKit✔️✔️:     res += '\n';
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut res = String::new();
    if let Some(name) = mol.properties().name() {
        if !name.is_empty() {
            res.push_str("COMPND    ");
            res.push_str(name);
            res.push('\n');
        }
    }

    // RDKit✔️✔️: unsigned int atm_count = 0;
    // RDKit✔️✔️: unsigned int ter_count = 0;
    // RDKit✔️✔️: unsigned int conect_count = 0;
    let mut atm_count: u32 = 0;
    let mut ter_count: u32 = 0;
    let mut conect_count: u32 = 0;

    // RDKit✔️✔️: const Conformer *conf;
    // RDKit✔️✔️: if (confId < 0 && rwmol.getNumConformers() > 1) {
    // RDKit✔️✔️:   int count = rwmol.getNumConformers();
    // RDKit✔️✔️:   for (confId = 0; confId < count; confId++) {
    // RDKit✔️✔️:     conf = &(rwmol.getConformer(confId));
    // RDKit✔️✔️:     ss << "MODEL     ";
    // RDKit✔️✔️:     ss << std::setw(4) << (confId + 1);
    // RDKit✔️✔️:     ss << "\n";
    // RDKit✔️✔️:     res += MolToPDBBody(rwmol, conf, flavor, atm_count, ter_count, conect_count);
    // RDKit✔️✔️:     res += "ENDMDL\n";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   if (confId < 0 && rwmol.getNumConformers() == 0) {
    // RDKit✔️✔️:     conf = nullptr;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     conf = &(rwmol.getConformer(confId));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let conformers = mol.conformers_3d();
    if conf_id < 0 && conformers.len() > 1 {
        // Write MODEL/ENDMDL for each conformer.
        for (idx, conf) in conformers.iter().enumerate() {
            res.push_str(&format!("MODEL     {:>4}\n", idx + 1));
            res.push_str(&mol_to_pdb_body(
                mol,
                Some(conf),
                flavor,
                &mut atm_count,
                &mut ter_count,
                &mut conect_count,
            ));
            res.push_str("ENDMDL\n");
        }
    } else {
        let conf = if conf_id < 0 {
            // No valid conformer selected; try first if available
            conformers.first()
        } else {
            conformers.get(conf_id as usize)
        };
        res.push_str(&mol_to_pdb_body(
            mol,
            conf,
            flavor,
            &mut atm_count,
            &mut ter_count,
            &mut conect_count,
        ));
    }

    // RDKit✔️✔️: if (flavor & 16) {
    // RDKit✔️✔️:   ss << "MASTER        0    0    0    0    0    0    0    0";
    // RDKit✔️✔️:   ss << std::setw(5) << atm_count;
    // RDKit✔️✔️:   ss << std::setw(5) << ter_count;
    // RDKit✔️✔️:   ss << std::setw(5) << conect_count;
    // RDKit✔️✔️:   ss << "    0\n";
    // RDKit✔️✔️:   res += ss.str();
    // RDKit✔️✔️: }
    if (flavor & 16) != 0 {
        res.push_str(&format!(
            "MASTER        0    0    0    0    0    0    0    0{:>5}{:>5}{:>5}    0\n",
            atm_count, ter_count, conect_count
        ));
    }

    // RDKit✔️✔️: res += "END\n";
    res.push_str("END\n");

    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;

    #[test]
    fn shared_periodic_table_delegation_covers_pdb_symbols_and_fallback() {
        for atomic_number in 0..=118 {
            assert_eq!(
                super::element_symbol(atomic_number),
                crate::rdkit_element_symbol(atomic_number).unwrap()
            );
        }
        assert_eq!(super::element_symbol(200), "?");
    }

    #[test]
    fn test_get_default_atom_number() {
        let mut elem: BTreeMap<u8, u32> = BTreeMap::new();
        // First carbon: should become "1 "
        let n1 = get_default_atom_number(6, &mut elem);
        assert_eq!(n1, "1 ");
        // Second carbon: "2 "
        let n2 = get_default_atom_number(6, &mut elem);
        assert_eq!(n2, "2 ");
    }

    #[test]
    fn test_mol_to_pdb_block_empty() {
        let mol = Molecule::new();
        let result = mol_to_pdb_block(&mol, -1, 0);
        assert!(result.ends_with("END\n"));
        assert_eq!(result, "END\n");
    }

    #[test]
    fn test_mol_to_pdb_block_simple() {
        // Use a molecule with an actual bond so default flavor emits CONECT.
        let mol = Molecule::from_smiles("CO").unwrap();
        let result = mol_to_pdb_block(&mol, -1, 0);
        assert!(result.contains("HETATM"));
        assert!(result.contains("END\n"));
        assert!(result.contains("CONECT"));
    }

    #[test]
    fn test_mol_to_pdb_block_with_name() {
        let mol = Molecule::from_smiles("CCO").unwrap().with_name("Ethanol");
        let result = mol_to_pdb_block(&mol, -1, 0);
        assert!(result.contains("COMPND    Ethanol"));
    }

    #[test]
    fn test_mol_to_pdb_block_no_conect() {
        let mol = Molecule::from_smiles("CO").unwrap();
        // flavor & 2 = don't write CONECT
        let result = mol_to_pdb_block(&mol, -1, 2);
        assert!(!result.contains("CONECT"));
    }

    #[test]
    fn test_mol_to_pdb_block_master() {
        let mol = Molecule::from_smiles("O").unwrap();
        // flavor & 16 = write MASTER record
        let result = mol_to_pdb_block(&mol, -1, 16);
        assert!(result.contains("MASTER"));
    }
}
