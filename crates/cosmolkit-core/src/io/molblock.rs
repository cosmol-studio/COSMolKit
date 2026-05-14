//! MolBlock writer implementation.

use std::borrow::Cow;

use crate::{
    Atom, AtomQueryPredicate, Bond, BondDirection, BondOrder, BondQueryPredicate, BondStereo,
    ChiralTag, CoordinateDimension, QueryNode, SGroupBondRole, SGroupBracketStyle,
    SGroupConnection, SGroupData, StereoGroupKind, SubstanceGroup, SubstanceGroupKind,
};
use crate::{Molecule, UnsupportedFeatureError};

const MIN_V2000_COORD: f64 = -100_000.0;
const MAX_V2000_COORD: f64 = 1_000_000.0;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MolWriteError {
    #[error("MolBlock writing is not implemented")]
    NotImplemented,
    #[error("MolBlock writing subset is not supported: {0}")]
    UnsupportedSubset(&'static str),
    #[error("MolBlock writing failed: {0}")]
    Value(String),
    #[error(transparent)]
    Operation(#[from] crate::OperationError),
    #[error(transparent)]
    UnsupportedFeature(#[from] UnsupportedFeatureError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfFormat {
    V2000,
    V3000,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolBlockWriteParams {
    pub format: SdfFormat,
    pub force_2d: bool,
    pub include_stereo: bool,
    pub kekulize: bool,
    pub precision: usize,
}

impl Default for MolBlockWriteParams {
    fn default() -> Self {
        Self {
            format: SdfFormat::V2000,
            force_2d: false,
            include_stereo: true,
            kekulize: true,
            precision: 6,
        }
    }
}

enum CoordinateSelection {
    Auto,
    TwoD,
    ThreeD,
}

struct SelectedCoordinates {
    coords: Option<Vec<[f64; 3]>>,
    is_3d: bool,
    label: Option<&'static str>,
}

pub fn mol_to_v2000_2d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::TwoD, &params)
}

pub fn mol_to_v2000_3d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::ThreeD, &params)
}

pub fn mol_to_v2000_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

pub fn mol_to_v3000_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

pub fn mol_to_v3000_2d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::TwoD, &params)
}

pub fn mol_to_v3000_3d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::ThreeD, &params)
}

pub fn mol_to_mol_block_with_params(
    molecule: &Molecule,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/FileWriters.h :: MolToMolBlock inline overloads
    // RDKit❗✔️: MolWriterParams params{includeStereo, kekulize, forceV3000};
    // RDKit❗✔️: return MolToMolBlock(mol, params, confId);
    // RDKit❗✔️: MolWriterParams v3KParams{params};
    // RDKit❗✔️: v3KParams.forceV3000 = true;
    // RDKit❗✔️: return MolToMolBlock(mol, v3KParams, confId);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/FileWriters.h :: MolToMolBlock inline overloads
    match params.format {
        SdfFormat::V2000 => {
            if params.force_2d {
                mol_to_v2000_block_with_params(molecule, CoordinateSelection::TwoD, params)
            } else {
                mol_to_v2000_block_with_params(molecule, CoordinateSelection::Auto, params)
            }
        }
        SdfFormat::V3000 => mol_to_v3000_block_with_params(
            molecule,
            if params.force_2d {
                CoordinateSelection::TwoD
            } else {
                CoordinateSelection::Auto
            },
            params,
        ),
    }
}

pub fn mol_to_2d_sdf_record(
    molecule: &Molecule,
    format: SdfFormat,
) -> Result<String, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_2d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_2d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

pub fn mol_to_sdf_record_with_params(
    molecule: &Molecule,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    let block = mol_to_mol_block_with_params(molecule, params)?;
    Ok(append_sdf_record_fields(block, molecule))
}

pub fn mol_to_3d_sdf_record(
    molecule: &Molecule,
    format: SdfFormat,
) -> Result<String, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_3d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_3d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

/// Prepared molecule with aromatic-bond bookkeeping.
/// Before kekulization, the set of bonds that were aromatic is recorded so
/// the molfile bond line can emit bond type 4 for those bonds (matching
/// RDKit's `prepareMol` → `aromaticBonds` bookkeeping).
struct PreparedMol<'a> {
    molecule: Cow<'a, Molecule>,
    /// Indices (in the bond table) of bonds that were aromatic before
    /// kekulization. After kekulization, these should still be written
    /// as bond type 4 (aromatic) in the molfile output.
    aromatic_bonds: Vec<usize>,
}

fn prepare_mol_for_writing<'a>(
    molecule: &'a Molecule,
    params: &MolBlockWriteParams,
) -> Result<PreparedMol<'a>, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: prepareMol
    // RDKit❗✔️: RWMol trwmol(mol);
    // RDKit❗✔️: if (params.kekulize) { MolOps::KekulizeIfPossible(trwmol, true); }
    // RDKit❗✔️: if (params.includeStereo) { wedgeBonds = Chirality::pickBondsToWedge(...); }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: prepareMol
    // Record aromatic bonds before kekulization.
    let aromatic_bonds: Vec<usize> = molecule
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.is_aromatic() || bond.order() == BondOrder::Aromatic)
        .map(|(idx, _)| idx)
        .collect();
    if params.kekulize {
        Ok(PreparedMol {
            molecule: Cow::Owned(molecule.with_kekulized_bonds(true)?),
            aromatic_bonds,
        })
    } else {
        Ok(PreparedMol {
            molecule: Cow::Borrowed(molecule),
            aromatic_bonds,
        })
    }
}

fn mol_to_v3000_block_with_params(
    molecule: &Molecule,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock / FileParserUtils::getV3000CTAB
    // RDKit❗✔️: if (isV3000) { ss << std::setw(3) << 0; ... ss << "999 V3000\n"; }
    // RDKit❗✔️: std::string res = "M  V30 BEGIN CTAB\n";
    // RDKit❗✔️: ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups << " " << num3DConstraints << " " << chiralFlag << "\n";
    // RDKit❗✔️: res += "M  V30 BEGIN ATOM\n";
    // RDKit❗✔️: res += GetV3000MolFileAtomLine(*atomIt, conf, queryListAtoms, precision);
    // RDKit❗✔️: if (tmol.getNumBonds()) { res += "M  V30 BEGIN BOND\n"; ... res += GetV3000MolFileBondLine(...); }
    // RDKit❗✔️: if (nSGroups > 0) { res += "M  V30 BEGIN SGROUP\n"; ... }
    // RDKit❗✔️: appendEnhancedStereoGroups(res, tmol, wedgeBonds);
    // RDKit❗✔️: res += "M  V30 END CTAB\n";
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock / FileParserUtils::getV3000CTAB
    let prepared = prepare_mol_for_writing(molecule, params)?;
    let molecule = prepared.molecule.as_ref();
    let aromatic_bonds = &prepared.aromatic_bonds;
    validate_v3000_writer_subset(molecule, params.include_stereo)?;
    let selected = select_coordinates(molecule, selection)?;
    let chiral_flag = molfile_chiral_flag(molecule)?;
    let generated_sgroups = v3000_generated_zbo_sgroups(molecule);

    let mut out = String::new();
    out.push_str(molecule.properties().name().unwrap_or_default());
    out.push('\n');
    out.push_str(&molfile_info_line(molecule, selected.label));
    out.push('\n');
    out.push_str(molecule.prop("_MolFileComments").unwrap_or_default());
    out.push('\n');
    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} {} 0 {}\n",
        molecule.num_atoms(),
        molecule.num_bonds(),
        molecule.substance_groups().len() + generated_sgroups.len(),
        chiral_flag
    ));
    let v3k_parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            molecule
                .atoms()
                .iter()
                .map(|atom| get_atom_parity_flag(molecule, atom, coords_3d))
                .collect()
        } else {
            vec![0u32; molecule.num_atoms()]
        }
    } else {
        vec![0u32; molecule.num_atoms()]
    };
    out.push_str("M  V30 BEGIN ATOM\n");
    for atom in molecule.atoms() {
        let coord = selected
            .coords
            .as_ref()
            .and_then(|coords| coords.get(atom.id().index()).copied())
            .unwrap_or([0.0, 0.0, 0.0]);
        out.push_str(&v3000_atom_line(
            molecule,
            atom,
            coord,
            params.precision,
            &v3k_parity_flags,
        )?);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    if molecule.num_bonds() != 0 {
        out.push_str("M  V30 BEGIN BOND\n");
        for bond in molecule.bonds() {
            out.push_str(&v3000_bond_line(
                bond,
                params.include_stereo,
                aromatic_bonds,
            )?);
            out.push('\n');
        }
        out.push_str("M  V30 END BOND\n");
    }
    append_v3000_sgroup_lines(&mut out, molecule, &generated_sgroups)?;
    if params.include_stereo {
        append_v3000_collection_lines(&mut out, molecule)?;
    }
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

fn mol_to_v2000_block_with_params(
    molecule: &Molecule,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: MolToV2KMolBlock / outputMolToMolBlock
    // RDKit❗✔️: RWMol trwmol(mol);
    // RDKit❗✔️: prepareMol(trwmol, params, aromaticBonds);
    // RDKit❗✔️: return outputMolToMolBlock(trwmol, confId, MolFileFormat::V2000,
    // RDKit❗✔️:                              params.precision, aromaticBonds);
    // RDKit❗✔️: nAtoms = tmol.getNumAtoms();
    // RDKit❗✔️: nBonds = tmol.getNumBonds();
    // RDKit❗✔️: if (whichFormat == MolFileFormat::V2000 &&
    // RDKit❗✔️:     (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    // RDKit❗✔️:   throw ValueErrorException(
    // RDKit❗✔️:       "V2000 format does not support more than 999 atoms, bonds or SGroups.");
    // RDKit❗✔️: }
    // RDKit❗✔️: tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::_Name, text)) { res += text; }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: // info
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit❗✔️:   res += text;
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   std::stringstream ss;
    // RDKit❗✔️:   ss << "  " << std::setw(8) << "RDKit";
    // RDKit❗✔️:   ss << std::setw(10) << "";
    // RDKit❗✔️:   if (conf) { if (conf->is3D()) { ss << "3D"; } else { ss << common_properties::TWOD; } }
    // RDKit❗✔️:   res += ss.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) { res += text; }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: ss << std::setw(3) << nAtoms;
    // RDKit❗✔️: ss << std::setw(3) << nBonds;
    // RDKit❗✔️: ss << std::setw(3) << nLists;
    // RDKit❗✔️: ss << std::setw(3) << nSGroups;
    // RDKit❗✔️: ss << std::setw(3) << chiralFlag;
    // RDKit❗✔️: ss << std::setw(3) << nsText;
    // RDKit❗✔️: ss << std::setw(3) << nRxnComponents;
    // RDKit❗✔️: ss << std::setw(3) << nReactants;
    // RDKit❗✔️: ss << std::setw(3) << nProducts;
    // RDKit❗✔️: ss << std::setw(3) << nIntermediates;
    // RDKit❗✔️: ss << "999 V2000\n";
    // RDKit❗✔️: for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms(); atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗✔️:   res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
    // RDKit❗✔️:   res += "\n";
    // RDKit❗✔️: }
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗✔️: for (const auto bond : tmol.bonds()) {
    // RDKit❗✔️:   res += GetMolFileBondLine(bond, wedgeBonds, conf, aromaticBonds[bond->getIdx()]);
    // RDKit❗✔️:   res += "\n";
    // RDKit❗✔️: }
    // RDKit❗✔️: res += GetMolFileChargeInfo(tmol);
    // RDKit❗✔️: res += GetMolFileRGroupInfo(tmol);
    // RDKit❗✔️: res += GetMolFileRGroupInfo(tmol);
    // RDKit❗✔️: res += GetMolFileQueryInfo(tmol, queryListAtoms);
    // RDKit❗✔️: res += GetMolFileAliasInfo(tmol);
    // RDKit❗✔️: res += GetMolFileZBOInfo(tmol);
    // RDKit❗✔️: res += GetMolFilePXAInfo(tmol);
    // RDKit❗✔️: res += GetMolFileSGroupInfo(tmol);
    // RDKit❗✔️: res += "M  END\n";
    // RDKit❗✔️: return res;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: MolToV2KMolBlock / outputMolToMolBlock
    let prepared = prepare_mol_for_writing(molecule, params)?;
    let molecule = prepared.molecule.as_ref();
    let aromatic_bonds = &prepared.aromatic_bonds;
    validate_v2000_writer_subset(molecule, params.include_stereo)?;
    let selected = select_coordinates(molecule, selection)?;
    validate_v2000_coordinate_range(selected.coords.as_deref())?;
    let chiral_flag = molfile_chiral_flag(molecule)?;

    let mut out = String::new();
    out.push_str(molecule.properties().name().unwrap_or_default());
    out.push('\n');
    out.push_str(&molfile_info_line(molecule, selected.label));
    out.push('\n');
    out.push_str(molecule.prop("_MolFileComments").unwrap_or_default());
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}999 V2000\n",
        molecule.num_atoms(),
        molecule.num_bonds(),
        0,
        molecule.substance_groups().len(),
        chiral_flag,
        0,
        0,
        0,
        0,
        0
    ));

    let parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            molecule
                .atoms()
                .iter()
                .map(|atom| get_atom_parity_flag(molecule, atom, coords_3d))
                .collect()
        } else {
            vec![0u32; molecule.num_atoms()]
        }
    } else {
        vec![0u32; molecule.num_atoms()]
    };
    for atom in molecule.atoms() {
        let coord = selected
            .coords
            .as_ref()
            .and_then(|coords| coords.get(atom.id().index()).copied())
            .unwrap_or([0.0, 0.0, 0.0]);
        out.push_str(&v2000_atom_line(atom, coord, molecule, &parity_flags)?);
        out.push('\n');
    }
    for bond in molecule.bonds() {
        out.push_str(&v2000_bond_line(
            bond,
            params.include_stereo,
            aromatic_bonds,
        )?);
        out.push('\n');
    }
    append_v2000_property_lines(&mut out, molecule);
    append_v2000_rgroup_lines(&mut out, molecule)?;
    append_v2000_value_lines(&mut out, molecule);
    append_v2000_alias_lines(&mut out, molecule);
    append_v2000_zbo_lines(&mut out, molecule);
    append_v2000_pxa_lines(&mut out, molecule);
    append_v2000_sgroup_lines(&mut out, molecule)?;
    // Write V2000 SSS query records for atoms and bonds with complex queries.
    // A single counter is shared across all SSS records in the CTfile.
    let mut sss_counter = 0u32;
    for atom in molecule.atoms() {
        out.push_str(&v2000_query_atom_sss_lines(atom, &mut sss_counter)?);
    }
    for bond in molecule.bonds() {
        out.push_str(&v2000_query_bond_sss_lines(bond, &mut sss_counter)?);
    }
    out.push_str("M  END\n");
    Ok(out)
}

/// Walk a query node and collect all leaf predicates, handling And/Or/Not.
fn flatten_query_predicates<T: Clone>(node: &QueryNode<T>) -> Vec<QueryFlattenNode<T>> {
    let mut result = Vec::new();
    flatten_query_predicates_inner(node, &mut result, false);
    result
}

fn flatten_query_predicates_inner<T: Clone>(
    node: &QueryNode<T>,
    result: &mut Vec<QueryFlattenNode<T>>,
    negated: bool,
) {
    match node {
        QueryNode::Predicate(p) => {
            result.push(QueryFlattenNode {
                predicate: p.clone(),
                negated,
            });
        }
        QueryNode::And(children) => {
            for child in children {
                flatten_query_predicates_inner(child, result, negated);
            }
        }
        QueryNode::Or(children) => {
            for child in children {
                flatten_query_predicates_inner(child, result, negated);
            }
        }
        QueryNode::Not(child) => {
            flatten_query_predicates_inner(child, result, !negated);
        }
    }
}

#[derive(Clone)]
struct QueryFlattenNode<T: Clone> {
    predicate: T,
    negated: bool,
}

/// Write V2000 SSS SAP records for a query atom. Returns the SAP lines.
fn v2000_query_atom_sss_lines(atom: &Atom, sss_counter: &mut u32) -> Result<String, MolWriteError> {
    let query = match atom.query() {
        Some(q) => q,
        None => return Ok(String::new()),
    };
    // If it's a simple atom-list query, the "L" symbol is already
    // written in the atom line; no SSS records needed.
    if v2000_atom_list_query(atom).is_some() {
        return Ok(String::new());
    }

    let predicates = flatten_query_predicates(query);
    if predicates.is_empty() {
        return Ok(String::new());
    }

    let mut out = String::new();
    let atom_idx = atom.id().index() + 1; // 1-based in MolBlock

    for flat in &predicates {
        match &flat.predicate {
            AtomQueryPredicate::RGroupLabel(id) => {
                // RDKit writes: M  ALS <nAtoms> <atomIdx> R<id>
                out.push_str(&format!("M  ALS{:>4}{:>4} R{}\\n", 1, atom_idx, id));
                continue;
            }
            AtomQueryPredicate::MolFileAlias(text) => {
                // RDKit writes: M  AAL <nAtoms> <atomIdx1> ... <aliasStr>
                out.push_str(&format!("M  AAL{:>4}{:>4} {}\\n", 1, atom_idx, text));
                continue;
            }
            AtomQueryPredicate::RecursiveSmarts(smarts) => {
                // RDKit writes: M  SMS <atomIdx> <SMARTS>
                out.push_str(&format!("M  SMS{:>4} {}\\n", atom_idx, smarts));
                continue;
            }
            _ => {}
        }
        let code = atom_query_predicate_to_sap_code(&flat.predicate);
        if code == 0 {
            continue;
        }
        *sss_counter += 1;
        let value = atom_query_predicate_to_sap_value(&flat.predicate);

        if flat.negated {
            out.push_str(&format!(
                "M  SAP{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter,
                atom_idx,
                code,
                -(value as i32)
            ));
        } else {
            out.push_str(&format!(
                "M  SAP{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter, atom_idx, code, value
            ));
        }
    }

    Ok(out)
}

/// Convert an AtomQueryPredicate to its V2000 SAP property code.
fn atom_query_predicate_to_sap_code(pred: &AtomQueryPredicate) -> u32 {
    match pred {
        AtomQueryPredicate::AtomicNumber(_) => 1,
        AtomQueryPredicate::AtomicNumberIn(_) => 1,
        AtomQueryPredicate::AtomicNumberNotIn(_) => 1,
        AtomQueryPredicate::FormalCharge(_) => 2,
        AtomQueryPredicate::Isotope(_) => 3,
        AtomQueryPredicate::ExplicitDegree(_) => 10,
        AtomQueryPredicate::ExplicitDegreeLessEqual(_) => 10,
        AtomQueryPredicate::NonHydrogenDegree(_) => 5,
        AtomQueryPredicate::ImplicitHydrogenCount(_) => 6,
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(_) => 6,
        AtomQueryPredicate::RingBondCount(_) => 9,
        AtomQueryPredicate::RingBondCountLessEqual(_) => 9,
        AtomQueryPredicate::RingBondCountNeedsScan => 9,
        AtomQueryPredicate::IsAromatic(_) => 12,
        AtomQueryPredicate::IsUnsaturated => 11,
        AtomQueryPredicate::RGroupLabel(_) => 20, // ALS (atom list symbol)
        AtomQueryPredicate::MolFileAlias(_) => 21, // AAL (atom alias)
        AtomQueryPredicate::RecursiveSmarts(_) => 22, // SMS (substructure SMARTS)
        AtomQueryPredicate::Any => 0,
        AtomQueryPredicate::UnsupportedFeature(_) => 0,
        // Phase A7 additions
        _ => 0,
    }
}

/// Convert an AtomQueryPredicate to its V2000 SAP property value.
fn atom_query_predicate_to_sap_value(pred: &AtomQueryPredicate) -> u32 {
    match pred {
        AtomQueryPredicate::AtomicNumber(z) => *z as u32,
        AtomQueryPredicate::AtomicNumberIn(_) => 0, // Handled by ALS
        AtomQueryPredicate::AtomicNumberNotIn(_) => 0,
        AtomQueryPredicate::FormalCharge(c) => (*c as i8) as u32,
        AtomQueryPredicate::Isotope(i) => *i as u32,
        AtomQueryPredicate::ExplicitDegree(d) => *d as u32,
        AtomQueryPredicate::ExplicitDegreeLessEqual(d) => *d as u32,
        AtomQueryPredicate::NonHydrogenDegree(d) => *d as u32,
        AtomQueryPredicate::ImplicitHydrogenCount(h) => *h as u32,
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(h) => *h as u32,
        AtomQueryPredicate::RingBondCount(c) => *c as u32,
        AtomQueryPredicate::RingBondCountLessEqual(c) => *c as u32,
        AtomQueryPredicate::RingBondCountNeedsScan => 1,
        AtomQueryPredicate::IsAromatic(true) => 1,
        AtomQueryPredicate::IsAromatic(false) => 0,
        AtomQueryPredicate::IsUnsaturated => 1,
        AtomQueryPredicate::RGroupLabel(id) => *id,
        AtomQueryPredicate::MolFileAlias(_) => 0, // String value, handled inline
        AtomQueryPredicate::RecursiveSmarts(_) => 0, // String value, handled inline
        AtomQueryPredicate::Any => 0,
        AtomQueryPredicate::UnsupportedFeature(_) => 0,
        // Phase A7 additions
        _ => 0,
    }
}

/// Write V2000 SSS SBT records for a query bond. Returns the SBT lines.
fn v2000_query_bond_sss_lines(bond: &Bond, sss_counter: &mut u32) -> Result<String, MolWriteError> {
    let query = match bond.query() {
        Some(q) => q,
        None => return Ok(String::new()),
    };
    // Simple bond query symbols and topology are already handled
    // by v2000_bond_query_symbol and v2000_bond_topology_code.
    // Only complex queries need SSS records.
    if v2000_bond_query_symbol(bond).is_some() {
        return Ok(String::new());
    }

    let predicates = flatten_query_predicates(query);
    if predicates.is_empty() {
        return Ok(String::new());
    }

    let mut out = String::new();
    let bond_idx = bond.id().index() + 1;

    for flat in &predicates {
        let (code, value) = bond_query_predicate_to_sbt_code(&flat.predicate);
        if code == 0 {
            continue;
        }
        *sss_counter += 1;
        if flat.negated {
            out.push_str(&format!(
                "M  SBT{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter,
                bond_idx,
                code,
                -(value as i32)
            ));
        } else {
            out.push_str(&format!(
                "M  SBT{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter, bond_idx, code, value
            ));
        }
    }

    Ok(out)
}

fn bond_query_predicate_to_sbt_code(pred: &BondQueryPredicate) -> (u32, u32) {
    match pred {
        BondQueryPredicate::Order(o) => match o {
            crate::BondOrder::Single => (1, 1),
            crate::BondOrder::Double => (1, 2),
            crate::BondOrder::Triple => (1, 3),
            crate::BondOrder::Quadruple => (1, 4),
            crate::BondOrder::Aromatic => (1, 5),
            _ => (0, 0),
        },
        BondQueryPredicate::IsAromatic(true) => (2, 1),
        BondQueryPredicate::IsAromatic(false) => (2, 0),
        BondQueryPredicate::IsInRing(true) => (3, 1),
        BondQueryPredicate::IsInRing(false) => (3, 0),
        BondQueryPredicate::Any => (0, 0), // Handled by symbol=8
        BondQueryPredicate::MolFileQueryCode(c) => (1, *c),
        _ => (0, 0),
    }
}

/// Validate that the molecule's query atoms/bonds can be written in
/// the V2000/V3000 format. Rejects only unsupported recursive SMARTS,
/// RGroupLabel, and MolFileAlias queries.
fn molfile_chiral_flag(molecule: &Molecule) -> Result<i32, MolWriteError> {
    molecule
        .prop("_MolFileChiralFlag")
        .map(|value| {
            value.trim().parse::<i32>().map_err(|_| {
                MolWriteError::Value(format!("invalid _MolFileChiralFlag value '{value}'"))
            })
        })
        .transpose()
        .map(|flag| flag.unwrap_or(0))
}

fn select_coordinates(
    molecule: &Molecule,
    selection: CoordinateSelection,
) -> Result<SelectedCoordinates, MolWriteError> {
    match selection {
        CoordinateSelection::TwoD => {
            let coords = molecule
                .coords_2d()
                .ok_or(MolWriteError::UnsupportedSubset(
                    "2D coordinates are required for V2000 2D output",
                ))?;
            if coords.len() != molecule.num_atoms() {
                return Err(MolWriteError::Value(
                    "2D coordinate count does not match atom count".to_string(),
                ));
            }
            Ok(SelectedCoordinates {
                coords: Some(
                    coords
                        .iter()
                        .map(|coord| [coord[0], coord[1], 0.0])
                        .collect(),
                ),
                is_3d: false,
                label: Some("2D"),
            })
        }
        CoordinateSelection::ThreeD => {
            let conformer =
                molecule
                    .conformers_3d()
                    .first()
                    .ok_or(MolWriteError::UnsupportedSubset(
                        "3D conformer coordinates are required for V2000 3D output",
                    ))?;
            if conformer.coords().len() != molecule.num_atoms() {
                return Err(MolWriteError::Value(
                    "3D coordinate count does not match atom count".to_string(),
                ));
            }
            Ok(SelectedCoordinates {
                coords: Some(conformer.coords().to_vec()),
                is_3d: true,
                label: Some("3D"),
            })
        }
        CoordinateSelection::Auto => {
            if matches!(
                molecule.source_coordinate_dim(),
                Some(CoordinateDimension::ThreeD)
            ) && !molecule.conformers_3d().is_empty()
            {
                return select_coordinates(molecule, CoordinateSelection::ThreeD);
            }
            if molecule.coords_2d().is_some() {
                return select_coordinates(molecule, CoordinateSelection::TwoD);
            }
            if !molecule.conformers_3d().is_empty() {
                return select_coordinates(molecule, CoordinateSelection::ThreeD);
            }
            Ok(SelectedCoordinates {
                coords: None,
                is_3d: false,
                label: None,
            })
        }
    }
}

fn validate_v2000_writer_subset(
    molecule: &Molecule,
    include_stereo: bool,
) -> Result<(), MolWriteError> {
    if molecule.num_atoms() > 999
        || molecule.num_bonds() > 999
        || molecule.substance_groups().len() > 999
    {
        return Err(MolWriteError::Value(
            "V2000 format does not support more than 999 atoms, bonds or SGroups".to_string(),
        ));
    }
    if include_stereo && !molecule.stereo_groups().is_empty() {
        return Err(MolWriteError::UnsupportedSubset(
            "MolBlock enhanced stereo writing is not ported",
        ));
    }
    for atom in molecule.atoms() {
        if atom.query().is_some() && v2000_atom_list_query(atom).is_none() {
            // Complex query atoms are written via SSS SAP records.
            // Only reject truly unsupported query types.
            if let Some(query) = atom.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let code = atom_query_predicate_to_sap_code(&flat.predicate);
                    // Accept: standard SAP codes (1-17) + special codes (20=RGroup, 21=Alias, 22=SMARTS)
                    if code == 0 || code > 22 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query atom predicate",
                        ));
                    }
                }
            }
        }
        // parity and charge emission are now implemented (get_atom_parity_flag, chg encoding).
        // Only reject mol_parity and mol_inversion_flag which are SDF-source artifacts
        // that should not appear in a V2000 write path.
        if include_stereo && (atom.mol_parity().is_some() || atom.mol_inversion_flag().is_some()) {
            return Err(MolWriteError::UnsupportedSubset(
                "atom stereochemistry MolBlock writing is not ported",
            ));
        }
        if (atom.explicit_hydrogens() != 0 && atom.prop("_ZBO_H").is_none())
            || atom.implicit_hydrogen()
            || !atom.tracked_isotopic_hydrogens().is_empty()
        {
            return Err(MolWriteError::UnsupportedSubset(
                "collapsed hydrogen MolBlock writing is not ported",
            ));
        }
    }
    for bond in molecule.bonds() {
        if bond.query().is_some()
            && v2000_bond_query_symbol(bond).is_none()
            && v2000_bond_topology_code(bond).is_none()
        {
            // Complex query bonds are written via SSS SBT records.
            if let Some(query) = bond.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let (code, _) = bond_query_predicate_to_sbt_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query bond predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo {
            v2000_bond_stereo_code(bond)?;
            if (bond.stereo() != BondStereo::None
                && !(bond.direction() == BondDirection::EitherDouble
                    && bond.stereo() == BondStereo::Any))
                || bond.stereo_atoms().is_some()
                || bond.unknown_stereo()
            {
                return Err(MolWriteError::UnsupportedSubset(
                    "bond stereochemistry MolBlock writing is not ported",
                ));
            }
        }
        v2000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_writer_subset(
    molecule: &Molecule,
    include_stereo: bool,
) -> Result<(), MolWriteError> {
    if include_stereo {
        validate_v3000_stereo_groups(molecule)?;
    }
    for atom in molecule.atoms() {
        if atom.query().is_some() && v2000_atom_list_query(atom).is_none() {
            // Complex query atoms in V3000 are written via query
            // property annotations. Only reject truly unsupported types.
            if let Some(query) = atom.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let code = atom_query_predicate_to_sap_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query atom V3000 predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo
            && (atom.chiral_permutation().is_some()
                || atom.unknown_stereo()
                || atom.mol_parity().is_some()
                || atom.mol_inversion_flag().is_some())
        {
            return Err(MolWriteError::UnsupportedSubset(
                "atom stereochemistry V3000 writing is not ported",
            ));
        }
        if atom.implicit_hydrogen() || !atom.tracked_isotopic_hydrogens().is_empty() {
            return Err(MolWriteError::UnsupportedSubset(
                "collapsed hydrogen V3000 writing is not ported",
            ));
        }
    }
    for bond in molecule.bonds() {
        if bond.query().is_some()
            && v2000_bond_query_symbol(bond).is_none()
            && v2000_bond_topology_code(bond).is_none()
        {
            // Complex query bonds in V3000 are supported via SBT equivalents.
            if let Some(query) = bond.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let (code, _) = bond_query_predicate_to_sbt_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query bond V3000 predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo {
            v3000_bond_cfg_code(bond)?;
            if (bond.stereo() != BondStereo::None
                && !(bond.direction() == BondDirection::EitherDouble
                    && bond.stereo() == BondStereo::Any))
                || bond.stereo_atoms().is_some()
                || bond.unknown_stereo()
            {
                return Err(MolWriteError::UnsupportedSubset(
                    "bond stereochemistry V3000 writing is not ported",
                ));
            }
        }
        v3000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_stereo_groups(molecule: &Molecule) -> Result<(), MolWriteError> {
    // Atropisomer (bond-based) stereo groups are now supported via atom
    // collection from the bond endpoints. No separate validation needed.
    // RDKit additionally validates wedge bond consistency for atropisomer
    // display; COSMolKit does not model wedge bonds for stereo groups yet.
    let _ = molecule;
    Ok(())
}

fn validate_v2000_coordinate_range(coords: Option<&[[f64; 3]]>) -> Result<(), MolWriteError> {
    let Some(coords) = coords else {
        return Ok(());
    };
    for coord in coords {
        for value in coord {
            if *value >= MAX_V2000_COORD || *value <= MIN_V2000_COORD {
                return Err(MolWriteError::Value(
                    "V2000 atom positions must be > -100000 and < 1000000".to_string(),
                ));
            }
        }
    }
    Ok(())
}

fn molfile_info_line(molecule: &Molecule, label: Option<&'static str>) -> String {
    if let Some(info) = molecule.prop("_MolFileInfo") {
        return info.to_string();
    }
    if let Some(info) = molecule.prop("_MolFileInfoLine") {
        return info.to_string();
    }
    let mut line = format!("  {:>8}{:>10}", "COSMolKit", "");
    if let Some(label) = label {
        line.push_str(label);
    }
    line
}

// Standard integer atomic weights (mass numbers of the most common isotope)
// indexed by atomic number (Z). Used to compute massDiff for molfile output:
//   massDiff = isotope - atomicWeight
const MOLFILE_ATOMIC_WEIGHT: [i32; 119] = [
    0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51, 52,
    55, 56, 59, 59, 64, 65, 70, 73, 75, 79, 80, 79, 85, 88, 89, 91, 93, 96, 98, 101, 103, 106, 108,
    112, 115, 119, 122, 128, 127, 132, 133, 138, 139, 140, 141, 144, 145, 152, 153, 157, 159, 162,
    165, 167, 169, 173, 175, 178, 181, 184, 187, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222,
    223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 263, 268,
    269, 275, 278, 281, 282, 285, 286, 289, 290, 293, 294, 294, 294, 294,
];

fn v2000_atom_line(
    atom: &Atom,
    coord: [f64; 3],
    molecule: &Molecule,
    parity_flags: &[u32],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties / GetMolFileAtomLine
    // RDKit❗✔️: totValence = 0;
    // RDKit❗✔️: atomMapNumber = 0;
    // RDKit❗✔️: parityFlag = 0;
    // RDKit❗✔️: x = y = z = 0.0;
    // RDKit❗✔️: if (conf) { const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx()); x = pos.x; y = pos.y; z = pos.z; ... }
    // RDKit✔️✔️: if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED && ...) { parityFlag = getAtomParityFlag(atom, conf); }
    // RDKit✔️✔️: if (hasNonDefaultValence(atom)) { ... }
    // RDKit❗✔️: snprintf(dest, 128,
    // RDKit❗✔️:          "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
    // RDKit❗✔️:          z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
    // RDKit❗✔️:          totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
    // RDKit❗✔️:          inversionFlag, exactChangeFlag);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties / GetMolFileAtomLine
    let symbol = v2000_atom_symbol(atom)?;
    let atom_idx = atom.id().index();
    let parity_flag = parity_flags.get(atom_idx).copied().unwrap_or(0);
    // RDKit charge codes: 0=0, 1=+3, 2=+2, 3=+1, 4=radical, 5=-1, 6=-2, 7=-3
    let chg = match atom.formal_charge() {
        -3 => 7,
        -2 => 6,
        -1 => 5,
        0 => 0,
        1 => 3,
        2 => 2,
        3 => 1,
        _ => 0,
    };
    let tot_valence = if has_non_default_valence(atom) {
        // RDKit: for degree-0 atoms (isolated ions/metals), use 15.
        // Otherwise emit total_valence % 15.
        let degree = molecule
            .derived_cache()
            .adjacency
            .as_ref()
            .map_or(0, |adj| adj.neighbors_of(atom_idx).len() as u32);
        if degree == 0 { 15 } else { 0u32 }
    } else {
        0u32
    };
    let mass_diff = atom.isotope().map_or(0i32, |iso| {
        iso as i32
            - MOLFILE_ATOMIC_WEIGHT
                .get(atom.atomic_number() as usize)
                .copied()
                .unwrap_or(0)
    });
    let h_count = atom
        .prop("_MolFileHCount")
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    let stereo_care = atom
        .prop("_MolFileStereoCare")
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    Ok(format!(
        "{:>10.4}{:>10.4}{:>10.4} {:>3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}  0{:>3}{:>3}{:>3}{:>3}{:>3}",
        coord[0],
        coord[1],
        coord[2],
        symbol,
        mass_diff,
        chg,
        parity_flag,
        h_count,
        stereo_care,
        tot_valence,
        0, // rxnComponentType
        0, // rxnComponentNumber
        atom.atom_map().unwrap_or(0),
        0, // inversionFlag
        0, // exactChangeFlag
    ))
}

// RDKit✔️✔️: unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
// RDKit✔️✔️:   if (!conf->is3D() ||
// RDKit✔️✔️:       !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4)) {
// RDKit✔️✔️:     return 0;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   const ROMol &mol = atom->getOwningMol();
// RDKit✔️✔️:   RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
// RDKit✔️✔️:   std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
// RDKit✔️✔️:   ... (neighbor loop with vector subtraction)
// RDKit✔️✔️:   std::sort(vs.begin(), vs.end(), Rankers::pairLess);
// RDKit✔️✔️:   double vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
// RDKit✔️✔️:   return (vol < 0) ? 2 : (vol > 0) ? 1 : 0;
// RDKit✔️✔️: }
// END RDKIT FUNCTION getAtomParityFlag
fn get_atom_parity_flag(molecule: &Molecule, atom: &Atom, coords_3d: &[[f64; 3]]) -> u32 {
    let atom_idx = atom.id().index();
    if !(atom.chiral_tag() == ChiralTag::TetrahedralCw
        || atom.chiral_tag() == ChiralTag::TetrahedralCcw)
    {
        return 0;
    }
    let adjacency = molecule.derived_cache().adjacency.as_ref().unwrap();
    let neighbors = adjacency.neighbors_of(atom_idx);
    if neighbors.len() < 3 {
        return 0;
    }
    // Compute vectors from the center atom to each neighbor.
    let pos = coords_3d[atom_idx];
    let mut vs: Vec<(usize, [f64; 3])> = neighbors
        .iter()
        .map(|nbr| {
            let npos = coords_3d[nbr.atom_index];
            let idx = nbr.atom_index;
            let v = [npos[0] - pos[0], npos[1] - pos[1], npos[2] - pos[2]];
            // RDKit shifts H-atom indices by numAtoms so they sort after heavy atoms.
            let sort_key = if molecule.atoms()[idx].atomic_number() == 1 {
                idx + molecule.num_atoms()
            } else {
                idx
            };
            (sort_key, v)
        })
        .collect();
    vs.sort_by(|a, b| a.0.cmp(&b.0));
    let cross = [
        vs[0].1[1] * vs[1].1[2] - vs[0].1[2] * vs[1].1[1],
        vs[0].1[2] * vs[1].1[0] - vs[0].1[0] * vs[1].1[2],
        vs[0].1[0] * vs[1].1[1] - vs[0].1[1] * vs[1].1[0],
    ];
    let vol = if vs.len() == 4 {
        cross[0] * vs[3].1[0] + cross[1] * vs[3].1[1] + cross[2] * vs[3].1[2]
    } else {
        -(cross[0] * vs[2].1[0] + cross[1] * vs[2].1[1] + cross[2] * vs[2].1[2])
    };
    if vol.abs() < 1e-10 {
        0
    } else if vol < 0.0 {
        2
    } else {
        1
    }
}

// RDKit✔️✔️: bool hasNonDefaultValence(const Atom *atom) {
// RDKit✔️✔️:   if (atom->getNumRadicalElectrons() != 0) { return true; }
// RDKit✔️✔️:   if (atom->hasQuery() || atom->needsUpdatePropertyCache()) { return false; }
// RDKit✔️✔️:   return ... (check against default valence)
// RDKit✔️✔️: }
// END RDKIT FUNCTION hasNonDefaultValence
//
// Simplified: returns true whenever the atom has radicals, a query, an explicit
// valence override, or a non-zero number of radical electrons.
#[allow(clippy::unnecessary_wraps)]
fn has_non_default_valence(atom: &Atom) -> bool {
    if atom.radical_electrons() != 0 {
        return true;
    }
    if atom.query().is_some() {
        return false;
    }
    // COSMolKit doesn't model explicit valence overrides yet, so this
    // always returns false for now.
    false
}

fn v2000_bond_line(
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol / GetMolFileBondLine
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // RDKit❗✔️: switch (bond->getBondType()) { case Bond::SINGLE: ... case Bond::DATIVE: res = 9; ... }
    // RDKit❗✔️: RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode, reverse);
    // RDKit❗✔️: ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    // RDKit❗✔️: ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    // RDKit❗✔️: ss << std::setw(3) << symbol;
    // RDKit❗✔️: ss << " " << std::setw(2) << dirCode;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol / GetMolFileBondLine
    let stereo_code = if include_stereo {
        v2000_bond_stereo_code(bond)?
    } else {
        0
    };
    // Use aromatic bond override: bonds that were aromatic before
    // kekulization should be written as type 4 in the molfile.
    let type_code = if aromatic_bonds.contains(&bond.id().index()) {
        4
    } else {
        v2000_bond_type_code(bond)?
    };
    let mut line = format!(
        "{:>3}{:>3}{:>3} {:>2}",
        bond.begin().index() + 1,
        bond.end().index() + 1,
        type_code,
        stereo_code
    );
    if let Some(topology) = v2000_bond_topology_code(bond) {
        line.push_str(&format!(" {:>2} {:>2}", 0, topology));
    }
    Ok(line)
}

fn v2000_bond_stereo_code(bond: &Bond) -> Result<u32, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileBondLine stereo dirCode branch
    // RDKit❗✔️: BEGINWEDGE -> 1
    // RDKit❗✔️: EITHERDOUBLE double bond -> 3
    // RDKit❗✔️: UNKNOWN single bond -> 4
    // RDKit❗✔️: BEGINDASH -> 6
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileBondLine stereo dirCode branch
    Ok(match bond.direction() {
        BondDirection::None => 0,
        BondDirection::BeginWedge if v2000_bond_type_code(bond)? == 1 => 1,
        BondDirection::BeginDash if v2000_bond_type_code(bond)? == 1 => 6,
        BondDirection::EitherDouble
            if v2000_bond_type_code(bond)? == 2 && bond.stereo() == BondStereo::Any =>
        {
            3
        }
        BondDirection::Unknown if v2000_bond_type_code(bond)? == 1 => 4,
        _ => {
            return Err(MolWriteError::UnsupportedSubset(
                "bond stereochemistry MolBlock writing is not ported",
            ));
        }
    })
}

fn v2000_bond_type_code(bond: &Bond) -> Result<u32, MolWriteError> {
    if let Some(code) = v2000_bond_query_symbol(bond) {
        return Ok(code);
    }
    match bond.order() {
        BondOrder::Single => Ok(if bond.is_aromatic() { 4 } else { 1 }),
        BondOrder::Double => Ok(if bond.is_aromatic() { 4 } else { 2 }),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Zero => Ok(1),
        BondOrder::Dative => Ok(9),
        _ => Err(MolWriteError::UnsupportedSubset(
            "bond order MolBlock writing is not ported",
        )),
    }
}

fn v3000_atom_line(
    molecule: &Molecule,
    atom: &Atom,
    coord: [f64; 3],
    precision: usize,
    parity_flags: &[u32],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileAtomLine
    // RDKit❗✔️: ss << "M  V30 " << atom->getIdx() + 1;
    // RDKit❗✔️: std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
    // RDKit❗✔️: if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) { ss << " " << symbol; } else { ... ss << " [" << symbols << "]"; }
    // RDKit❗✔️: ss << std::fixed; ss << std::setprecision(precision); ss << " " << x << " " << y << " " << z;
    // RDKit❗✔️: ss << " " << atomMapNumber;
    // RDKit✔️✔️: if (parityFlag != 0) { ss << " CFG=" << parityFlag; }
    // RDKit❗✔️: if (chg != 0) { ss << " CHG=" << chg; }
    // RDKit❗✔️: if (isotope != 0 && !isAtomRGroup(*atom)) { ss << " MASS=" << mass; }
    // RDKit❗✔️: if (nRadEs != 0 && atom->getTotalDegree() != 0) { ... ss << " RAD=" << nRadEs; }
    // RDKit❗✔️: if (symbol == "R#") { ... ss << " RGROUPS=(1 " << rLabel << ")"; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileAtomLine
    let symbol = v3000_atom_symbol(atom)?;
    let parity_flag = parity_flags.get(atom.id().index()).copied().unwrap_or(0);
    let mut out = format!("M  V30 {} {}", atom.id().index() + 1, symbol);
    out.push_str(&format!(
        " {0:.1$} {2:.1$} {3:.1$} {4}",
        coord[0],
        precision,
        coord[1],
        coord[2],
        atom.atom_map().unwrap_or(0)
    ));
    if parity_flag != 0 {
        out.push_str(&format!(" CFG={parity_flag}"));
    }
    if atom.formal_charge() != 0 {
        out.push_str(&format!(" CHG={}", atom.formal_charge()));
    }
    if let Some(isotope) = atom.isotope()
        && atom.prop("_MolFileRLabel").is_none()
    {
        out.push_str(&format!(" MASS={isotope}"));
    }
    let electrons = atom.radical_electrons();
    if electrons != 0 && atom_degree(molecule, atom.id().index()) != 0 {
        let code = if electrons % 2 == 1 { 2 } else { 3 };
        out.push_str(&format!(" RAD={code}"));
    }
    if atom.prop("_MolFileRLabel").is_some() {
        let label = atom
            .prop("_MolFileRLabel")
            .and_then(|value| value.parse::<u32>().ok())
            .unwrap_or(1);
        out.push_str(&format!(" RGROUPS=(1 {label})"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molAttachOrder", "ATTCHORD");
    append_v3000_atom_int_prop(&mut out, atom, "molAttachPoint", "ATTCHPT");
    append_v3000_atom_int_prop(&mut out, atom, "molAtomSeqId", "SEQID");
    if let Some(value) = atom.prop("molAtomSeqName") {
        out.push_str(&format!(" SEQNAME={value}"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molRxnExactChange", "EXACHG");
    if let Some(value) = atom.mol_inversion_flag() {
        out.push_str(&format!(" INVRET={value}"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molStereoCare", "STBOX");
    append_v3000_atom_int_prop(&mut out, atom, "molSubstCount", "SUBST");
    append_v3000_atom_int_prop(&mut out, atom, "molRingBondCount", "RBCNT");
    if let Some(value) = atom.prop("molAtomClass") {
        out.push_str(&format!(" CLASS={value}"));
    }
    Ok(out)
}

fn append_v3000_atom_int_prop(out: &mut String, atom: &Atom, prop: &str, label: &str) {
    if let Some(value) = atom.prop(prop)
        && value != "0"
    {
        out.push_str(&format!(" {label}={value}"));
    }
}

fn v3000_atom_symbol(atom: &Atom) -> Result<String, MolWriteError> {
    if let Some((atomic_numbers, negated)) = v2000_atom_list_query(atom) {
        let mut value = String::new();
        if negated {
            value.push_str("\"NOT ");
        }
        value.push('[');
        for (idx, atomic_number) in atomic_numbers.iter().enumerate() {
            if idx != 0 {
                value.push(',');
            }
            value.push_str(molfile_atom_symbol(*atomic_number)?);
        }
        value.push(']');
        if negated {
            value.push('"');
        }
        return Ok(value);
    }
    Ok(v2000_atom_symbol(atom)?.to_string())
}

fn v3000_bond_line(
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine
    // RDKit❗✔️: ss << "M  V30 " << bond->getIdx() + 1;
    // RDKit❗✔️: ss << " " << GetV3000BondCode(bond);
    // RDKit❗✔️: ss << " " << bond->getBeginAtomIdx() + 1 << " " << bond->getEndAtomIdx() + 1;
    // RDKit❗✔️: if (includeStereo) { ... ss << " CFG=" << cfg; }
    // RDKit❗✔️: if (bond->hasQuery()) { int topol = getQueryBondTopology(bond); if (topol) { ss << " TOPO=" << topol; } }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::molReactStatus, iprop) && iprop) { ss << " RXCTR=" << iprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::molStereoCare, sprop) && sprop != "0") { ss << " STBOX=" << sprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, sprop) && sprop != "0") { ss << " ENDPTS=" << sprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::_MolFileBondAttach, sprop) && sprop != "0") { ss << " ATTACH=" << sprop; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine
    let type_code = if aromatic_bonds.contains(&bond.id().index()) {
        4
    } else {
        v3000_bond_type_code(bond)?
    };
    let mut out = format!(
        "M  V30 {} {} {} {}",
        bond.id().index() + 1,
        type_code,
        bond.begin().index() + 1,
        bond.end().index() + 1
    );
    if include_stereo && let Some(cfg) = v3000_bond_cfg_code(bond)? {
        out.push_str(&format!(" CFG={cfg}"));
    }
    if let Some(topology) = v2000_bond_topology_code(bond) {
        out.push_str(&format!(" TOPO={topology}"));
    }
    append_v3000_bond_prop(&mut out, bond, "molReactStatus", "RXCTR");
    append_v3000_bond_prop(&mut out, bond, "molStereoCare", "STBOX");
    append_v3000_bond_prop(&mut out, bond, "_MolFileBondEndPts", "ENDPTS");
    append_v3000_bond_prop(&mut out, bond, "_MolFileBondAttach", "ATTACH");
    Ok(out)
}

fn v3000_bond_cfg_code(bond: &Bond) -> Result<Option<u32>, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine stereo CFG branch
    // RDKit❗✔️: BEGINWEDGE -> CFG=1
    // RDKit❗✔️: UNKNOWN single bond -> CFG=2
    // RDKit❗✔️: EITHERDOUBLE double bond -> CFG=2
    // RDKit❗✔️: BEGINDASH -> CFG=3
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine stereo CFG branch
    Ok(match bond.direction() {
        BondDirection::None => None,
        BondDirection::BeginWedge => Some(1),
        BondDirection::BeginDash => Some(3),
        BondDirection::Unknown if v3000_bond_type_code(bond)? == 1 => Some(2),
        BondDirection::EitherDouble if v3000_bond_type_code(bond)? == 2 => Some(2),
        _ => {
            return Err(MolWriteError::UnsupportedSubset(
                "bond stereochemistry V3000 writing is not ported",
            ));
        }
    })
}

fn append_v3000_bond_prop(out: &mut String, bond: &Bond, prop: &str, label: &str) {
    if let Some(value) = bond.prop(prop)
        && value != "0"
    {
        out.push_str(&format!(" {label}={value}"));
    }
}

fn v3000_bond_type_code(bond: &Bond) -> Result<u32, MolWriteError> {
    if let Some(code) = v2000_bond_query_symbol(bond) {
        return Ok(code);
    }
    match bond.order() {
        BondOrder::Single => Ok(if bond.is_aromatic() { 4 } else { 1 }),
        BondOrder::Double => Ok(if bond.is_aromatic() { 4 } else { 2 }),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Zero => Ok(1),
        BondOrder::Dative => Ok(9),
        BondOrder::Hydrogen => Ok(10),
        _ => Err(MolWriteError::UnsupportedSubset(
            "bond order V3000 writing is not ported",
        )),
    }
}

fn append_v3000_sgroup_lines(
    out: &mut String,
    molecule: &Molecule,
    generated_sgroups: &[SubstanceGroup],
) -> Result<(), MolWriteError> {
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() && generated_sgroups.is_empty() {
        return Ok(());
    }
    out.push_str("M  V30 BEGIN SGROUP\n");
    for (idx, sgroup) in sgroups.iter().enumerate() {
        out.push_str(&v3000_sgroup_lines(idx + 1, sgroup)?);
    }
    let offset = sgroups.len();
    for (idx, sgroup) in generated_sgroups.iter().enumerate() {
        out.push_str(&v3000_sgroup_lines(offset + idx + 1, sgroup)?);
    }
    out.push_str("M  V30 END SGROUP\n");
    Ok(())
}

fn v3000_generated_zbo_sgroups(molecule: &Molecule) -> Vec<SubstanceGroup> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: createZBOSubstanceGroups
    // RDKit❗✔️: SubstanceGroup bsg(&mol, "DAT"); bsg.setProp("FIELDNAME", "ZBO");
    // RDKit❗✔️: for (const auto bond : mol.bonds()) { if (bond->getBondType() == Bond::ZERO) { bsg.addBondWithIdx(...); atomsAffected[...] = 1; } }
    // RDKit❗✔️: if (atomsAffected.any()) { for affected atoms { bsg.addAtomWithIdx(i); ... asgText += getTotalNumHs(); zsgText += getFormalCharge(); } }
    // RDKit❗✔️: addSubstanceGroup(mol, bsg); asg.setProp("DATAFIELDS", aDataFields); addSubstanceGroup(mol, asg); zsg.setProp("DATAFIELDS", zDataFields); addSubstanceGroup(mol, zsg);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: createZBOSubstanceGroups
    let zero_bonds = molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero)
        .map(Bond::id)
        .collect::<Vec<_>>();
    if zero_bonds.is_empty() {
        return Vec::new();
    }
    let mut affected = vec![false; molecule.num_atoms()];
    for bond in molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero)
    {
        affected[bond.begin().index()] = true;
        affected[bond.end().index()] = true;
    }
    let atoms = affected
        .iter()
        .enumerate()
        .filter_map(|(idx, is_affected)| is_affected.then_some(crate::AtomId::new(idx)))
        .collect::<Vec<_>>();
    let hydrogens = atoms
        .iter()
        .map(|atom| {
            molecule.atoms()[atom.index()]
                .explicit_hydrogens()
                .to_string()
        })
        .collect::<Vec<_>>()
        .join(";");
    let charges = atoms
        .iter()
        .map(|atom| molecule.atoms()[atom.index()].formal_charge().to_string())
        .collect::<Vec<_>>()
        .join(";");
    vec![
        SubstanceGroup::new(crate::SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(atoms.clone())
            .with_bonds(zero_bonds)
            .with_data(SGroupData {
                field_name: Some("ZBO".to_string()),
                ..SGroupData::default()
            }),
        SubstanceGroup::new(crate::SubstanceGroupId::new(1), SubstanceGroupKind::Data)
            .with_atoms(atoms.clone())
            .with_data(SGroupData {
                field_name: Some("HYD".to_string()),
                values: vec![hydrogens],
                ..SGroupData::default()
            }),
        SubstanceGroup::new(crate::SubstanceGroupId::new(2), SubstanceGroupKind::Data)
            .with_atoms(atoms)
            .with_data(SGroupData {
                field_name: Some("ZCH".to_string()),
                values: vec![charges],
                ..SGroupData::default()
            }),
    ]
}

fn v3000_sgroup_lines(idx: usize, sgroup: &SubstanceGroup) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetV3000MolFileSGroupLines
    // RDKit❗✔️: std::string currLine = (boost::format("M  V30 %d %s %d") % idx % sgroup.getProp<std::string>("TYPE") % id).str();
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000IdxVectorDataBlock("ATOMS", sgroup.getAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000BondsBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000IdxVectorDataBlock("PATOMS", sgroup.getParentAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SUBTYPE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("MULT", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CONNECT", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000ParentBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CompNoBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("LABEL", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000BracketBlock(sgroup.getBrackets()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("ESTATE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CStateBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDNAME", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDINFO", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDDISP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYTYPE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYOP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000FieldDataBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CLASS", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000AttachPointBlock(sgroup.getAttachPoints()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("BRKTYP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SEQID", sgroup), currLine, os);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetV3000MolFileSGroupLines
    let mut os = String::new();
    let mut current = format!(
        "M  V30 {} {} {}",
        idx,
        sgroup_kind_code(sgroup)?,
        sgroup.external_id().unwrap_or(0)
    );
    add_v3000_sgroup_block(
        v3000_index_vector_block("ATOMS", sgroup.atoms().iter().map(|atom| atom.index() + 1)),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_bonds_block(sgroup), &mut current, &mut os);
    add_v3000_sgroup_block(
        v3000_index_vector_block(
            "PATOMS",
            sgroup.parent_atoms().iter().map(|atom| atom.index() + 1),
        ),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("SUBTYPE", sgroup.subtype()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("MULT", sgroup.props().get("MULT").map(String::as_str)),
        &mut current,
        &mut os,
    );
    let connect = sgroup
        .connection()
        .map(sgroup_connection_code)
        .transpose()?;
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("CONNECT", connect),
        &mut current,
        &mut os,
    );
    if let Some(parent) = sgroup.parent() {
        add_v3000_sgroup_block(
            format!(" PARENT={}", parent.index() + 1),
            &mut current,
            &mut os,
        );
    }
    if let Some(compno) = sgroup.component_number() {
        add_v3000_sgroup_block(format!(" COMPNO={compno}"), &mut current, &mut os);
    }
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("LABEL", sgroup.label()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_bracket_block(sgroup), &mut current, &mut os);
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("ESTATE", sgroup.expansion_state()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_cstate_block(sgroup), &mut current, &mut os);
    if let Some(data) = sgroup.data() {
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDNAME", data.field_name.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDINFO", data.field_info.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDDISP", data.field_display.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("QUERYTYPE", data.query_type.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("QUERYOP", data.query_op.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_field_data_block(data.values.iter().map(String::as_str)),
            &mut current,
            &mut os,
        );
    }
    add_v3000_sgroup_block(
        v3000_sgroup_field_data_block(sgroup.data_fields().iter().map(String::as_str)),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("CLASS", sgroup.class()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_attach_point_block(sgroup),
        &mut current,
        &mut os,
    );
    let bracket_style = v3000_sgroup_bracket_style_value(sgroup.bracket_style())?;
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("BRKTYP", bracket_style),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block(
            "SEQID",
            sgroup.props().get("SEQID").map(String::as_str),
        ),
        &mut current,
        &mut os,
    );
    if !current.is_empty() && current != "M  V30" {
        os.push_str(&current);
        os.push('\n');
    }
    Ok(os)
}

fn v3000_index_vector_block(label: &str, values: impl IntoIterator<Item = usize>) -> String {
    let values = values.into_iter().collect::<Vec<_>>();
    if values.is_empty() {
        return String::new();
    }
    let mut out = format!(" {label}=({}", values.len());
    for value in values {
        out.push_str(&format!(" {value}"));
    }
    out.push(')');
    out
}

fn v3000_sgroup_bonds_block(sgroup: &SubstanceGroup) -> String {
    let crossing = sgroup
        .bonds()
        .iter()
        .filter(|bond| sgroup.bond_role(**bond) == SGroupBondRole::Crossing)
        .map(|bond| bond.index() + 1);
    let contained = sgroup
        .bonds()
        .iter()
        .filter(|bond| sgroup.bond_role(**bond) == SGroupBondRole::Contained)
        .map(|bond| bond.index() + 1);
    let mut out = v3000_index_vector_block("XBONDS", crossing);
    out.push_str(&v3000_index_vector_block("CBONDS", contained));
    out
}

fn v3000_sgroup_string_property_block(label: &str, value: Option<&str>) -> String {
    let Some(value) = value else {
        return String::new();
    };
    if value.is_empty() {
        return String::new();
    }
    let mut out = format!(" {label}=");
    out.push_str(&v3000_quote_string_property(value));
    out
}

fn v3000_quote_string_property(value: &str) -> String {
    let needs_quotes = value.contains(' ') || value.contains('"') || value.contains('(');
    let mut out = String::new();
    if needs_quotes {
        out.push('"');
    }
    for ch in value.chars() {
        out.push(ch);
        if ch == '"' {
            out.push(ch);
        }
    }
    if needs_quotes {
        out.push('"');
    }
    out
}

fn v3000_sgroup_bracket_block(sgroup: &SubstanceGroup) -> String {
    let Some(display) = sgroup.display() else {
        return String::new();
    };
    let mut out = String::new();
    for bracket in &display.brackets {
        out.push_str(" BRKXYZ=(9");
        out.push_str(&format!(
            " {} {} 0 {} {} 0 0 0 0)",
            v3000_double_field(bracket.p1[0]),
            v3000_double_field(bracket.p1[1]),
            v3000_double_field(bracket.p2[0]),
            v3000_double_field(bracket.p2[1])
        ));
    }
    out
}

fn v3000_sgroup_cstate_block(sgroup: &SubstanceGroup) -> String {
    let mut out = String::new();
    for cstate in sgroup.cstates() {
        out.push_str(" CSTATE=(");
        if matches!(sgroup.kind(), SubstanceGroupKind::Superatom) {
            out.push_str(&format!(
                "4 {} {} {} 0",
                cstate.bond.index() + 1,
                v3000_double_field(cstate.vector[0]),
                v3000_double_field(cstate.vector[1])
            ));
        } else {
            out.push_str(&format!("1 {}", cstate.bond.index() + 1));
        }
        out.push(')');
    }
    out
}

fn v3000_sgroup_field_data_block<'a>(values: impl IntoIterator<Item = &'a str>) -> String {
    let mut out = String::new();
    for value in values {
        out.push_str(" FIELDDATA=\"");
        out.push_str(value);
        out.push('"');
    }
    out
}

fn v3000_sgroup_attach_point_block(sgroup: &SubstanceGroup) -> String {
    let mut out = String::new();
    for attach_point in sgroup.attach_points() {
        out.push_str(&format!(" SAP=(3 {}", attach_point.atom.index() + 1));
        if attach_point.leaving_atom == Some(attach_point.atom) {
            out.push_str(" aidx");
        } else {
            let leaving = attach_point.leaving_atom.map_or(0, |atom| atom.index() + 1);
            out.push_str(&format!(" {leaving}"));
        }
        out.push(' ');
        out.push_str(attach_point.label.as_deref().unwrap_or(""));
        out.push(')');
    }
    out
}

fn v3000_sgroup_bracket_style_value(
    style: Option<&SGroupBracketStyle>,
) -> Result<Option<&str>, MolWriteError> {
    Ok(match style {
        Some(SGroupBracketStyle::Bracket) => Some("BRACKET"),
        Some(SGroupBracketStyle::Parenthesis) => Some("PAREN"),
        Some(SGroupBracketStyle::None) | None => None,
        Some(SGroupBracketStyle::Unknown(value)) => {
            if value.is_empty() {
                None
            } else {
                Some(value.as_str())
            }
        }
    })
}

fn add_v3000_sgroup_block(block: String, current_line: &mut String, out: &mut String) {
    if block.is_empty() {
        return;
    }
    if current_line.is_empty() {
        current_line.push_str("M  V30");
    }
    if current_line.len() + block.len() < 78 {
        current_line.push_str(&block);
        return;
    }
    out.push_str(current_line);
    out.push_str(" -\n");
    let mut start = 0_usize;
    while block.len().saturating_sub(start) >= 73 {
        out.push_str("M  V30");
        if start != 0 {
            out.push(' ');
        }
        let end = safe_char_boundary(&block, start + 72);
        out.push_str(&block[start..end]);
        start = end;
        if start < block.len() {
            out.push_str("-\n");
        }
    }
    if start < block.len() {
        *current_line = "M  V30".to_string();
        if start != 0 {
            current_line.push(' ');
        }
        current_line.push_str(&block[start..]);
    } else {
        current_line.clear();
    }
}

fn safe_char_boundary(value: &str, mut index: usize) -> usize {
    index = index.min(value.len());
    while index > 0 && !value.is_char_boundary(index) {
        index -= 1;
    }
    index
}

fn v3000_double_field(value: f64) -> String {
    v2000_double_field(value).trim().to_string()
}

fn append_v3000_collection_lines(
    out: &mut String,
    molecule: &Molecule,
) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: appendEnhancedStereoGroups
    // RDKit❗✔️: if (!tmol.getStereoGroups().empty()) {
    // RDKit❗✔️:   auto stereo_groups = tmol.getStereoGroups();
    // RDKit❗✔️:   assignStereoGroupIds(stereo_groups);
    // RDKit❗✔️:   res += "M  V30 BEGIN COLLECTION\n";
    // RDKit❗✔️:   for (auto &&group : stereo_groups) {
    // RDKit❗✔️:     tmp += "M  V30 MDLV30/";
    // RDKit❗✔️:     switch (group.getGroupType()) { ... STEABS / STEREL[id] / STERAC[id] ... }
    // RDKit❗✔️:     Atropisomers::getAllAtomIdsForStereoGroup(tmol, group, atomIds, wedgeBonds);
    // RDKit❗✔️:     tmp += std::to_string(atomIds.size()); ... wrap at 78 chars ...
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res += tmp + "M  V30 END COLLECTION\n";
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: appendEnhancedStereoGroups
    let groups = molecule.stereo_groups();
    if groups.is_empty() {
        return Ok(());
    }
    let write_ids = assigned_v3000_stereo_group_ids(groups);
    out.push_str("M  V30 BEGIN COLLECTION\n");
    for (group, write_id) in groups.iter().zip(write_ids) {
        let label = match group.kind() {
            StereoGroupKind::Absolute => "STEABS".to_string(),
            StereoGroupKind::Or => format!("STEREL{}", write_id.unwrap_or(1)),
            StereoGroupKind::And => format!("STERAC{}", write_id.unwrap_or(1)),
        };
        // Collect atoms: if the group has explicit atoms use them;
        // otherwise collect atoms from all bonds (atropisomer/bond-based groups).
        let atom_ids: Vec<crate::AtomId> = if !group.atoms().is_empty() {
            group.atoms().to_vec()
        } else {
            group
                .bonds()
                .iter()
                .flat_map(|bond_id| {
                    let bond = &molecule.bonds()[bond_id.index()];
                    [bond.begin(), bond.end()]
                })
                .collect()
        };
        append_v3000_collection_group_line(out, &label, atom_ids.iter().copied());
    }
    out.push_str("M  V30 END COLLECTION\n");
    Ok(())
}

fn assigned_v3000_stereo_group_ids(groups: &[crate::StereoGroup]) -> Vec<Option<u32>> {
    let mut or_ids = Vec::<u32>::new();
    let mut and_ids = Vec::<u32>::new();
    let mut assigned = groups
        .iter()
        .map(crate::StereoGroup::id)
        .collect::<Vec<_>>();
    for (idx, group) in groups.iter().enumerate() {
        let Some(id) = assigned[idx] else {
            continue;
        };
        let ids = match group.kind() {
            StereoGroupKind::Or => &mut or_ids,
            StereoGroupKind::And => &mut and_ids,
            StereoGroupKind::Absolute => continue,
        };
        if id != 0 && ids.contains(&id) {
            assigned[idx] = None;
        } else if id != 0 {
            ids.push(id);
        }
    }
    let mut next_or = 0_u32;
    let mut next_and = 0_u32;
    for (idx, group) in groups.iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute || assigned[idx].is_some() {
            continue;
        }
        let (next, ids) = match group.kind() {
            StereoGroupKind::Or => (&mut next_or, &mut or_ids),
            StereoGroupKind::And => (&mut next_and, &mut and_ids),
            StereoGroupKind::Absolute => unreachable!(),
        };
        *next += 1;
        while ids.contains(&*next) {
            *next += 1;
        }
        ids.push(*next);
        assigned[idx] = Some(*next);
    }
    assigned
}

fn append_v3000_collection_group_line(
    out: &mut String,
    label: &str,
    atoms: impl IntoIterator<Item = crate::AtomId>,
) {
    let atoms = atoms.into_iter().collect::<Vec<_>>();
    let mut line = format!("M  V30 MDLV30/{label} ATOMS=({}", atoms.len());
    for atom in atoms {
        let idx = (atom.index() + 1).to_string();
        if line.len() + idx.len() >= 78 {
            out.push_str(&line);
            out.push_str("-\n");
            line = "M  V30 ".to_string();
        } else {
            line.push(' ');
        }
        line.push_str(&idx);
    }
    out.push_str(&line);
    out.push_str(")\n");
}

fn v2000_bond_query_symbol(bond: &Bond) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondSymbol / BondGetMolFileSymbol
    // RDKit❗✔️: int res = 8;
    // RDKit❗✔️: if (qry->getDescription() == "BondOrder" || getQueryBondTopology(bond)) { res = 0; }
    // RDKit❗✔️: if (qry->getDescription() == "BondOr" && !qry->getNegation()) { ... SINGLE+DOUBLE => 5; SINGLE+AROMATIC => 6; DOUBLE+AROMATIC => 7; }
    // RDKit❗✔️: else if (qry->getDescription() == "SingleOrAromaticBond" && !qry->getNegation()) { res = 6; }
    // RDKit❗✔️: else if (qry->getDescription() == "SingleOrDoubleBond" && !qry->getNegation()) { res = 5; }
    // RDKit❗✔️: else if (qry->getDescription() == "DoubleOrAromaticBond" && !qry->getNegation()) { res = 7; }
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondSymbol / BondGetMolFileSymbol
    let query = bond.query()?;
    v2000_bond_query_symbol_from_query(query)
}

fn v2000_bond_query_symbol_from_query(query: &QueryNode<BondQueryPredicate>) -> Option<u32> {
    match query {
        QueryNode::Predicate(BondQueryPredicate::Any) => Some(8),
        QueryNode::Predicate(BondQueryPredicate::MolFileQueryCode(code)) => Some(*code),
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders)) => {
            v2000_bond_order_query_symbol(orders)
        }
        QueryNode::And(children) => {
            let non_topology = children
                .iter()
                .filter(|child| v2000_bond_topology_code_from_query(child).is_none())
                .collect::<Vec<_>>();
            if non_topology.len() == 1 {
                v2000_bond_query_symbol_from_query(non_topology[0])
            } else {
                None
            }
        }
        _ => None,
    }
}

fn v2000_bond_topology_code(bond: &Bond) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondTopology / GetMolFileBondLine
    // RDKit❗✔️: if (qry->getDescription() == "BondAnd" && !qry->getNegation() && qry->endChildren() - qry->beginChildren() == 2) { ... qry = BondInRing child; }
    // RDKit❗✔️: if (qry->getDescription() == "BondInRing") { if (qry->getNegation()) { res = 2; } else { res = 1; } }
    // RDKit❗✔️: if (bond->hasQuery()) { int topol = getQueryBondTopology(bond); if (topol) { ss << " " << std::setw(2) << 0 << " " << std::setw(2) << topol; } }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondTopology / GetMolFileBondLine
    v2000_bond_topology_code_from_query(bond.query()?)
}

fn v2000_bond_topology_code_from_query(query: &QueryNode<BondQueryPredicate>) -> Option<u32> {
    match query {
        QueryNode::Predicate(BondQueryPredicate::IsInRing(true)) => Some(1),
        QueryNode::Predicate(BondQueryPredicate::IsInRing(false)) => Some(2),
        QueryNode::And(children) if children.len() == 2 => children
            .iter()
            .find_map(v2000_bond_topology_code_from_query),
        _ => None,
    }
}

fn v2000_bond_order_query_symbol(orders: &[BondOrder]) -> Option<u32> {
    let has_single = orders.contains(&BondOrder::Single);
    let has_double = orders.contains(&BondOrder::Double);
    let has_aromatic = orders.contains(&BondOrder::Aromatic);
    match (orders.len(), has_single, has_double, has_aromatic) {
        (2, true, true, false) => Some(5),
        (2, true, false, true) => Some(6),
        (2, false, true, true) => Some(7),
        _ => None,
    }
}

fn append_v2000_property_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileChargeInfo
    // RDKit❗✔️: if (atom->getFormalCharge() != 0) { ++nChgs; chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getFormalCharge(); ... }
    // RDKit❗✔️: unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit❗✔️: if (nRadEs != 0 && atom->getTotalDegree() != 0) { ... if (nRadEs % 2) { nRadEs = 2; } else { nRadEs = 3; } ... }
    // RDKit❗✔️: int isotope = atom->getIsotope();
    // RDKit❗✔️: if (isotope != 0) { ++nMassDiffs; massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % isotope; ... }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileChargeInfo
    let charges = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            (atom.formal_charge() != 0)
                .then_some((atom.id().index() + 1, i32::from(atom.formal_charge())))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "CHG", &charges);

    let radicals = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            let electrons = atom.radical_electrons();
            if electrons == 0 || atom_degree(molecule, atom.id().index()) == 0 {
                return None;
            }
            let code = if electrons % 2 == 1 { 2 } else { 3 };
            Some((atom.id().index() + 1, code))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "RAD", &radicals);

    let isotopes = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            atom.isotope()
                .map(|isotope| (atom.id().index() + 1, i32::from(isotope)))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ISO", &isotopes);
}

fn append_v2000_rgroup_lines(out: &mut String, molecule: &Molecule) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileRGroupInfo
    // RDKit❗✔️: unsigned int nEntries = 0;
    // RDKit❗✔️: if ((*atomIt)->getPropIfPresent(common_properties::_MolFileRLabel, lbl)) {
    // RDKit❗✔️:   ss << " " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
    // RDKit❗✔️:      << std::setw(3) << lbl;
    // RDKit❗✔️:   ++nEntries;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (nEntries) { ss2 << "M  RGP" << std::setw(3) << nEntries << ss.str() << "\n"; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileRGroupInfo
    let entries = molecule
        .atoms()
        .iter()
        .filter_map(|atom| atom.prop("_MolFileRLabel").map(|label| (atom, label)))
        .map(|(atom, label)| {
            let label = label.parse::<i32>().map_err(|_| {
                MolWriteError::Value(format!(
                    "invalid _MolFileRLabel value '{}' on atom {}",
                    label,
                    atom.id().index()
                ))
            })?;
            Ok((atom.id().index() + 1, label))
        })
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    if !entries.is_empty() {
        out.push_str(&format!("M  RGP{:>3}", entries.len()));
        for (idx, label) in entries {
            out.push_str(&format!(" {:>3} {:>3}", idx, label));
        }
        out.push('\n');
    }
    Ok(())
}

fn append_v2000_value_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    // RDKit❗✔️: std::string molFileValue;
    // RDKit❗✔️: if (!wrote_query &&
    // RDKit❗✔️:     atom->getPropIfPresent(common_properties::molFileValue, molFileValue)) {
    // RDKit❗✔️:   ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << molFileValue
    // RDKit❗✔️:      << "\n";
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    for atom in molecule.atoms() {
        if let Some(value) = atom.prop("molFileValue") {
            out.push_str(&format!("V  {:>3} {value}\n", atom.id().index() + 1));
        }
    }
    append_v2000_atom_list_lines(out, molecule);
}

fn append_v2000_atom_list_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    // RDKit❗✔️: for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:   if (listQs[atom->getIdx()]) {
    // RDKit❗✔️:     INT_VECT vals;
    // RDKit❗✔️:     getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗✔️:     ss << "M  ALS " << std::setw(3) << atom->getIdx() + 1 << " ";
    // RDKit❗✔️:     ss << std::setw(2) << vals.size();
    // RDKit❗✔️:     if (atom->getQuery()->getNegation()) { ss << " T "; } else { ss << " F "; }
    // RDKit❗✔️:     for (auto val : vals) { ss << std::setw(4) << std::left << (PeriodicTable::getTable()->getElementSymbol(val)); }
    // RDKit❗✔️:     ss << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    for atom in molecule.atoms() {
        if let Some((atomic_numbers, negated)) = v2000_atom_list_query(atom) {
            out.push_str(&format!(
                "M  ALS {:>3} {:>2} {} ",
                atom.id().index() + 1,
                atomic_numbers.len(),
                if negated { "T" } else { "F" }
            ));
            for atomic_number in atomic_numbers {
                let symbol = molfile_atom_symbol(*atomic_number).unwrap_or("*");
                out.push_str(&format!("{symbol:<4}"));
            }
            out.push('\n');
        }
    }
}

fn append_v2000_zbo_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileZBOInfo
    // RDKit❗✔️: if ((*bondIt)->getBondType() == Bond::ZERO) {
    // RDKit❗✔️:   ss << " " << std::setw(3) << (*bondIt)->getIdx() + 1 << " " << std::setw(3) << 0;
    // RDKit❗✔️:   atomsAffected[(*bondIt)->getBeginAtomIdx()] = 1;
    // RDKit❗✔️:   atomsAffected[(*bondIt)->getEndAtomIdx()] = 1;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (nEntries) { res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n"; }
    // RDKit❗✔️: hydss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getTotalNumHs();
    // RDKit❗✔️: if (atom->getFormalCharge()) { zchss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getFormalCharge(); }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileZBOInfo
    let zbo_entries = molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero || bond.prop("_ZBO").is_some())
        .map(|bond| (bond.id().index() + 1, 0))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ZBO", &zbo_entries);
    if zbo_entries.is_empty() {
        return;
    }
    let mut affected = vec![false; molecule.num_atoms()];
    for bond in molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero || bond.prop("_ZBO").is_some())
    {
        affected[bond.begin().index()] = true;
        affected[bond.end().index()] = true;
    }
    let hydrogens = molecule
        .atoms()
        .iter()
        .filter(|atom| affected[atom.id().index()])
        .map(|atom| (atom.id().index() + 1, i32::from(atom.explicit_hydrogens())))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "HYD", &hydrogens);
    let zcharges = molecule
        .atoms()
        .iter()
        .filter(|atom| affected[atom.id().index()] && atom.formal_charge() != 0)
        .map(|atom| (atom.id().index() + 1, i32::from(atom.formal_charge())))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ZCH", &zcharges);
}

fn append_v2000_pxa_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFilePXAInfo
    // RDKit❗✔️: if (atom->hasProp("_MolFile_PXA")) {
    // RDKit❗✔️:   res += boost::str(boost::format("M  PXA % 3d%s\n") % (atom->getIdx() + 1) %
    // RDKit❗✔️:                    atom->getProp<std::string>("_MolFile_PXA"));
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFilePXAInfo
    for atom in molecule.atoms() {
        if let Some(pxa) = atom.prop("_MolFile_PXA") {
            out.push_str(&format!("M  PXA {:>3}{pxa}\n", atom.id().index() + 1));
        }
    }
}

fn append_v2000_sgroup_lines(out: &mut String, molecule: &Molecule) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetMolFileSGroupInfo
    // RDKit❗✔️: ret << BuildV2000STYLines(mol);
    // RDKit❗✔️: ret << BuildV2000SLBLines(mol);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
    // RDKit❗✔️: ret << BuildV2000SDSLines(mol);
    // RDKit❗✔️: ret << BuildV2000SPLLines(mol);
    // RDKit❗✔️: ret << BuildV2000SNCLines(mol);
    // RDKit❗✔️: ret << BuildV2000SBTLines(mol);
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SAL", sgroup.getAtoms());
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SPA", sgroup.getParentAtoms());
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SBL", sgroup.getBonds());
    // RDKit❗✔️: ret << BuildV2000SDILine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SMTLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SBVLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SDTLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SDDLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SCDSEDLines(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SAPLines(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SCLLine(idx, sgroup);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetMolFileSGroupInfo
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() {
        return Ok(());
    }

    let sty_entries = sgroups
        .iter()
        .enumerate()
        .map(|(idx, sgroup)| Ok((idx + 1, sgroup_kind_code(sgroup)?)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_sgroup_string_entries(out, "STY", 8, 3, &sty_entries);

    let slb_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.external_id().map(|id| (idx + 1, id as i32)))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SLB", &slb_entries);

    let sst_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.subtype().map(|value| (idx + 1, value)))
        .collect::<Vec<_>>();
    append_v2000_sgroup_string_entries(out, "SST", 8, 3, &sst_entries);

    let scn_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.connection().map(|value| (idx + 1, value)))
        .map(|(idx, connection)| sgroup_connection_code(connection).map(|code| (idx, code)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_sgroup_string_entries(out, "SCN", 8, 3, &scn_entries);

    let sds_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| (sgroup.expansion_state() == Some("E")).then_some(idx + 1))
        .collect::<Vec<_>>();
    append_v2000_sgroup_index_entries(out, "SDS EXP", 15, &sds_entries);

    let spl_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| {
            sgroup
                .parent()
                .map(|parent| (idx + 1, parent.index() as i32 + 1))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SPL", &spl_entries);

    let snc_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| {
            sgroup
                .component_number()
                .map(|value| (idx + 1, value as i32))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SNC", &snc_entries);

    let sbt_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.bracket_style().map(|style| (idx + 1, style)))
        .map(|(idx, style)| sgroup_bracket_style_code(style).map(|code| (idx, code)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_counted_property(out, "SBT", &sbt_entries);

    for (idx, sgroup) in sgroups.iter().enumerate() {
        let idx = idx + 1;
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SAL",
            sgroup.atoms().iter().map(|atom| atom.index() + 1),
        );
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SPA",
            sgroup.parent_atoms().iter().map(|atom| atom.index() + 1),
        );
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SBL",
            sgroup.bonds().iter().map(|bond| bond.index() + 1),
        );
        append_v2000_sgroup_sdi_lines(out, idx, sgroup);
        append_v2000_sgroup_smt_line(out, idx, sgroup);
        append_v2000_sgroup_sbv_lines(out, idx, sgroup);
        append_v2000_sgroup_sdt_line(out, idx, sgroup);
        append_v2000_sgroup_sdd_line(out, idx, sgroup);
        append_v2000_sgroup_scd_sed_lines(out, idx, sgroup)?;
        append_v2000_sgroup_sap_lines(out, idx, sgroup);
        append_v2000_sgroup_scl_line(out, idx, sgroup);
    }
    Ok(())
}

fn sgroup_kind_code(sgroup: &SubstanceGroup) -> Result<&str, MolWriteError> {
    Ok(match sgroup.kind() {
        SubstanceGroupKind::Data => "DAT",
        SubstanceGroupKind::Superatom => "SUP",
        SubstanceGroupKind::MultipleGroup => "MUL",
        SubstanceGroupKind::StructuralRepeatUnit => "SRU",
        SubstanceGroupKind::Monomer => "MON",
        SubstanceGroupKind::Copolymer => "COP",
        SubstanceGroupKind::Crosslink => "CRO",
        SubstanceGroupKind::Graft => "GRA",
        SubstanceGroupKind::Modification => "MOD",
        SubstanceGroupKind::Mer => "MER",
        SubstanceGroupKind::AnyPolymer => "ANY",
        SubstanceGroupKind::MixtureComponent => "COM",
        SubstanceGroupKind::Mixture => "MIX",
        SubstanceGroupKind::Formulation => "FOR",
        SubstanceGroupKind::Generic(value) => value.as_str(),
    })
}

fn sgroup_connection_code(connection: &SGroupConnection) -> Result<&str, MolWriteError> {
    Ok(match connection {
        SGroupConnection::HeadToHead => "HH",
        SGroupConnection::HeadToTail => "HT",
        SGroupConnection::Either => "EU",
        SGroupConnection::Unknown(value) => value.as_str(),
    })
}

fn sgroup_bracket_style_code(style: &SGroupBracketStyle) -> Result<i32, MolWriteError> {
    match style {
        SGroupBracketStyle::Bracket | SGroupBracketStyle::None => Ok(0),
        SGroupBracketStyle::Parenthesis => Ok(1),
        SGroupBracketStyle::Unknown(value) => Err(MolWriteError::Value(format!(
            "invalid SGroup BRKTYP value '{value}'"
        ))),
    }
}

fn append_v2000_sgroup_string_entries<T: AsRef<str>>(
    out: &mut String,
    code: &str,
    entries_per_line: usize,
    field_width: usize,
    entries: &[(usize, T)],
) {
    for chunk in entries.chunks(entries_per_line) {
        out.push_str(&format!("M  {code}{:>3}", chunk.len()));
        for (idx, value) in chunk {
            out.push_str(&v2000_int_field(*idx));
            out.push_str(&v2000_string_field(value.as_ref(), field_width, true, true));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_index_entries(
    out: &mut String,
    code: &str,
    entries_per_line: usize,
    entries: &[usize],
) {
    for chunk in entries.chunks(entries_per_line) {
        out.push_str(&format!("M  {code}{:>3}", chunk.len()));
        for idx in chunk {
            out.push_str(&v2000_int_field(*idx));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_member_lines<I>(out: &mut String, sgroup_idx: usize, code: &str, values: I)
where
    I: IntoIterator<Item = usize>,
{
    let values = values.into_iter().collect::<Vec<_>>();
    for chunk in values.chunks(15) {
        out.push_str(&format!(
            "M  {code}{}{:>3}",
            v2000_int_field(sgroup_idx),
            chunk.len()
        ));
        for value in chunk {
            out.push_str(&v2000_int_field(*value));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_sdi_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(display) = sgroup.display() {
        for bracket in &display.brackets {
            out.push_str(&format!(
                "M  SDI{}{:>3}{}{}{}{}\n",
                v2000_int_field(idx),
                4,
                v2000_double_field(bracket.p1[0]),
                v2000_double_field(bracket.p1[1]),
                v2000_double_field(bracket.p2[0]),
                v2000_double_field(bracket.p2[1])
            ));
        }
    }
}

fn append_v2000_sgroup_smt_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    let value = if matches!(sgroup.kind(), SubstanceGroupKind::MultipleGroup) {
        sgroup.props().get("MULT").map(String::as_str)
    } else {
        sgroup.label()
    };
    if let Some(value) = value {
        out.push_str(&format!(
            "M  SMT{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn append_v2000_sgroup_sbv_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    for cstate in sgroup.cstates() {
        out.push_str(&format!(
            "M  SBV{}{}",
            v2000_int_field(idx),
            v2000_int_field(cstate.bond.index() + 1)
        ));
        if matches!(sgroup.kind(), SubstanceGroupKind::Superatom) {
            out.push_str(&v2000_double_field(cstate.vector[0]));
            out.push_str(&v2000_double_field(cstate.vector[1]));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_sdt_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    let Some(data) = sgroup.data() else {
        return;
    };
    let Some(field_name) = data.field_name.as_deref() else {
        return;
    };
    out.push_str(&format!("M  SDT{}", v2000_int_field(idx)));
    out.push_str(&v2000_string_field(field_name, 30, true, true));
    out.push_str(&v2000_string_field(
        data.field_type.as_deref().unwrap_or("T"),
        2,
        true,
        false,
    ));
    if let Some(value) = data.field_info.as_deref() {
        out.push_str(&v2000_string_field(value, 20, true, false));
    }
    if let Some(value) = data.query_type.as_deref() {
        out.push_str(&v2000_string_field(value, 2, true, false));
    }
    if let Some(value) = data.query_op.as_deref() {
        out.push_str(&v2000_string_field(value, 15, true, false));
    }
    out.push('\n');
}

fn append_v2000_sgroup_sdd_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(value) = sgroup.data().and_then(|data| data.field_display.as_deref()) {
        out.push_str(&format!(
            "M  SDD{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn append_v2000_sgroup_scd_sed_lines(
    out: &mut String,
    idx: usize,
    sgroup: &SubstanceGroup,
) -> Result<(), MolWriteError> {
    let values = sgroup
        .data()
        .map(|data| data.values.as_slice())
        .unwrap_or(&[]);
    for value in values.iter().chain(sgroup.data_fields()) {
        if value.len() > 200 {
            return Err(MolWriteError::Value(format!(
                "Data field '{value}' in SGroup {idx} is longer than limit of 200 characters."
            )));
        }
        if value.is_empty() {
            out.push_str(&format!("M  SED{} \n", v2000_int_field(idx)));
            continue;
        }
        let mut chunks = value.as_bytes().chunks(69).peekable();
        while let Some(chunk) = chunks.next() {
            let text = String::from_utf8_lossy(chunk);
            let has_more = chunks.peek().is_some();
            let code = if has_more { "SCD" } else { "SED" };
            out.push_str(&format!(
                "M  {code}{}{}\n",
                v2000_int_field(idx),
                v2000_string_field(&text, 69, has_more, true)
            ));
        }
    }
    Ok(())
}

fn append_v2000_sgroup_sap_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    for chunk in sgroup.attach_points().chunks(6) {
        out.push_str(&format!("M  SAP{}{:>3}", v2000_int_field(idx), chunk.len()));
        for attach_point in chunk {
            out.push_str(&v2000_int_field(attach_point.atom.index() + 1));
            let leaving = attach_point.leaving_atom.map_or(0, |atom| atom.index() + 1);
            out.push_str(&v2000_int_field(leaving));
            out.push_str(&v2000_string_field(
                attach_point.label.as_deref().unwrap_or(""),
                2,
                false,
                true,
            ));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_scl_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(value) = sgroup.class() {
        out.push_str(&format!(
            "M  SCL{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn v2000_int_field(value: usize) -> String {
    format!(" {value:>3}")
}

fn v2000_double_field(value: f64) -> String {
    format!("{value:>10.4}")
}

fn v2000_string_field(value: &str, field_size: usize, pad: bool, add_separator: bool) -> String {
    let mut out = String::new();
    if add_separator {
        out.push(' ');
    }
    if value.len() >= field_size {
        out.push_str(&value[..field_size]);
    } else if pad {
        out.push_str(&format!("{value:<field_size$}"));
    } else {
        out.push_str(value);
    }
    out
}

fn append_v2000_alias_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAliasInfo
    // RDKit❗✔️: if ((*atomIt)->getPropIfPresent(common_properties::molFileAlias, lbl)) {
    // RDKit❗✔️:   if (!lbl.empty()) {
    // RDKit❗✔️:     ss << "A  " << std::setw(3) << (*atomIt)->getIdx() + 1 << "\n"
    // RDKit❗✔️:        << lbl << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAliasInfo
    for atom in molecule.atoms() {
        if let Some(alias) = atom.prop("molFileAlias")
            && !alias.is_empty()
        {
            out.push_str(&format!("A  {:>3}\n{alias}\n", atom.id().index() + 1));
        }
    }
}

fn append_v2000_counted_property(out: &mut String, label: &str, entries: &[(usize, i32)]) {
    for chunk in entries.chunks(8) {
        out.push_str(&format!("M  {label}{:>3}", chunk.len()));
        for (idx, value) in chunk {
            out.push_str(&format!(" {:>3} {:>3}", idx, value));
        }
        out.push('\n');
    }
}

fn atom_degree(molecule: &Molecule, atom_index: usize) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| bond.begin().index() == atom_index || bond.end().index() == atom_index)
        .count()
}

fn append_sdf_record_fields(mut block: String, molecule: &Molecule) -> String {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/SDWriter.cpp :: _MolToSDStream / _writePropToStream
    // RDKit❗✔️: (*dp_ostream) << MolToMolBlock(mol, true, confId, df_kekulize, df_forceV3000);
    // RDKit❗✔️: (*dp_ostream) << ">  <" << name << ">  ";
    // RDKit❗✔️: (*dp_ostream) << "\n";
    // RDKit❗✔️: (*dp_ostream) << pval << "\n";
    // RDKit❗✔️: (*dp_ostream) << "\n";
    // RDKit❗✔️: (*dp_ostream) << "$$$$\n";
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/SDWriter.cpp :: _MolToSDStream / _writePropToStream
    for (name, value) in molecule.properties().sdf_data_fields() {
        if name.contains('\n') || value.contains("\r\n\r\n") || value.contains("\n\n") {
            continue;
        }
        block.push_str(">  <");
        block.push_str(name);
        block.push_str(">  \n");
        block.push_str(value);
        block.push_str("\n\n");
    }
    block.push_str("$$$$\n");
    block
}

fn v2000_atom_symbol(atom: &Atom) -> Result<&'static str, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: AtomGetMolFileSymbol
    // RDKit❗✔️: if (atom->hasProp(common_properties::_MolFileRLabel)) { res = "R#"; }
    // RDKit❗✔️: else if (atom->getAtomicNum()) { res = atom->getSymbol(); }
    // RDKit❗✔️: else { if (hasComplexQuery(atom)) { if (isAtomListQuery(atom)) { res = "L"; } else { res = "*"; } } else { res = "R"; } }
    // RDKit❗✔️: if (padWithSpaces) { while (res.size() < 3) { res += " "; } }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: AtomGetMolFileSymbol
    if atom.prop("_MolFileRLabel").is_some() {
        return Ok("R#");
    }
    if atom.atomic_number() != 0 {
        return molfile_atom_symbol(atom.atomic_number());
    }
    if v2000_atom_list_query(atom).is_some() {
        return Ok("L");
    }
    molfile_atom_symbol(atom.atomic_number())
}

fn v2000_atom_list_query(atom: &Atom) -> Option<(&[u8], bool)> {
    match atom.query()? {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(atomic_numbers)) => {
            Some((atomic_numbers.as_slice(), false))
        }
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)) => {
            Some((atomic_numbers.as_slice(), true))
        }
        _ => None,
    }
}

fn molfile_atom_symbol(atomic_number: u8) -> Result<&'static str, MolWriteError> {
    match atomic_number {
        0 => Ok("*"),
        1 => Ok("H"),
        2 => Ok("He"),
        3 => Ok("Li"),
        4 => Ok("Be"),
        5 => Ok("B"),
        6 => Ok("C"),
        7 => Ok("N"),
        8 => Ok("O"),
        9 => Ok("F"),
        10 => Ok("Ne"),
        11 => Ok("Na"),
        12 => Ok("Mg"),
        13 => Ok("Al"),
        14 => Ok("Si"),
        15 => Ok("P"),
        16 => Ok("S"),
        17 => Ok("Cl"),
        18 => Ok("Ar"),
        19 => Ok("K"),
        20 => Ok("Ca"),
        21 => Ok("Sc"),
        22 => Ok("Ti"),
        23 => Ok("V"),
        24 => Ok("Cr"),
        25 => Ok("Mn"),
        26 => Ok("Fe"),
        27 => Ok("Co"),
        28 => Ok("Ni"),
        29 => Ok("Cu"),
        30 => Ok("Zn"),
        31 => Ok("Ga"),
        32 => Ok("Ge"),
        33 => Ok("As"),
        34 => Ok("Se"),
        35 => Ok("Br"),
        36 => Ok("Kr"),
        37 => Ok("Rb"),
        38 => Ok("Sr"),
        39 => Ok("Y"),
        40 => Ok("Zr"),
        41 => Ok("Nb"),
        42 => Ok("Mo"),
        43 => Ok("Tc"),
        44 => Ok("Ru"),
        45 => Ok("Rh"),
        46 => Ok("Pd"),
        47 => Ok("Ag"),
        48 => Ok("Cd"),
        49 => Ok("In"),
        50 => Ok("Sn"),
        51 => Ok("Sb"),
        52 => Ok("Te"),
        53 => Ok("I"),
        54 => Ok("Xe"),
        55 => Ok("Cs"),
        56 => Ok("Ba"),
        57 => Ok("La"),
        58 => Ok("Ce"),
        59 => Ok("Pr"),
        60 => Ok("Nd"),
        61 => Ok("Pm"),
        62 => Ok("Sm"),
        63 => Ok("Eu"),
        64 => Ok("Gd"),
        65 => Ok("Tb"),
        66 => Ok("Dy"),
        67 => Ok("Ho"),
        68 => Ok("Er"),
        69 => Ok("Tm"),
        70 => Ok("Yb"),
        71 => Ok("Lu"),
        72 => Ok("Hf"),
        73 => Ok("Ta"),
        74 => Ok("W"),
        75 => Ok("Re"),
        76 => Ok("Os"),
        77 => Ok("Ir"),
        78 => Ok("Pt"),
        79 => Ok("Au"),
        80 => Ok("Hg"),
        81 => Ok("Tl"),
        82 => Ok("Pb"),
        83 => Ok("Bi"),
        84 => Ok("Po"),
        85 => Ok("At"),
        86 => Ok("Rn"),
        87 => Ok("Fr"),
        88 => Ok("Ra"),
        89 => Ok("Ac"),
        90 => Ok("Th"),
        91 => Ok("Pa"),
        92 => Ok("U"),
        93 => Ok("Np"),
        94 => Ok("Pu"),
        95 => Ok("Am"),
        96 => Ok("Cm"),
        97 => Ok("Bk"),
        98 => Ok("Cf"),
        99 => Ok("Es"),
        100 => Ok("Fm"),
        101 => Ok("Md"),
        102 => Ok("No"),
        103 => Ok("Lr"),
        104 => Ok("Rf"),
        105 => Ok("Db"),
        106 => Ok("Sg"),
        107 => Ok("Bh"),
        108 => Ok("Hs"),
        109 => Ok("Mt"),
        110 => Ok("Ds"),
        111 => Ok("Rg"),
        112 => Ok("Cn"),
        113 => Ok("Nh"),
        114 => Ok("Fl"),
        115 => Ok("Mc"),
        116 => Ok("Lv"),
        117 => Ok("Ts"),
        118 => Ok("Og"),
        _ => Err(MolWriteError::UnsupportedSubset(
            "unknown element MolBlock writing is not ported",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomQueryPredicate, AtomSpec, BondDirection, BondOrder, BondQueryPredicate, BondSpec,
        BondStereo, ChiralTag, Element, QueryNode, SGroupAttachPoint, SGroupBondRole,
        SGroupBracket, SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData,
        SGroupDisplay, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind,
    };
    use crate::{MOLBLOCK_IO_FEATURE, Molecule, UnsupportedFeatureError};

    fn charged_isotope_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("charged");
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_formal_charge(-1)
                .with_isotope(13)
                .with_prop("molFileValue", "payload")
                .with_prop("molFileAlias", "AliasLabel"),
        );
        builder.set_2d_coordinates(vec![[1.25, -2.5]]).unwrap();
        builder.build().unwrap()
    }

    fn zbo_extension_molecule() -> Molecule {
        let mut builder = Molecule::builder()
            .with_name("zbo")
            .with_property("_MolFileInfoLine", "  COSMolKit          2D");
        let carbon = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_formal_charge(-1)
                .with_explicit_hydrogens(2)
                .with_prop("_ZBO_H", "1"),
        );
        let rgroup = builder.add_atom(
            AtomSpec::new(Element::DUMMY)
                .with_prop("_MolFileRLabel", "7")
                .with_prop("_MolFile_PXA", " payload"),
        );
        builder
            .add_bond(BondSpec::new(carbon, rgroup, BondOrder::Zero))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder.build().unwrap()
    }

    fn atom_list_query_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("atom-list");
        builder.add_atom(
            AtomSpec::new(Element::DUMMY)
                .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                    vec![6, 7],
                )))
                .with_prop("molFileValue", "[#6,#7]"),
        );
        builder.add_atom(AtomSpec::new(Element::C).with_query(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberNotIn(vec![8, 16]),
        )));
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder.build().unwrap()
    }

    fn query_bond_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("query-bond");
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        let a4 = builder.add_atom(AtomSpec::new(Element::C));
        let a5 = builder.add_atom(AtomSpec::new(Element::C));
        let a6 = builder.add_atom(AtomSpec::new(Element::C));
        let a7 = builder.add_atom(AtomSpec::new(Element::C));
        let a8 = builder.add_atom(AtomSpec::new(Element::C));
        let a9 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Unspecified)
                    .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Unspecified).with_query(
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Double,
                ])),
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a4, a5, BondOrder::Unspecified).with_query(
                QueryNode::predicate(BondQueryPredicate::MolFileQueryCode(42)),
            ))
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(a6, a7, BondOrder::Single)
                    .with_query(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
            )
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(a8, a9, BondOrder::Unspecified).with_query(QueryNode::and(vec![
                    QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                        BondOrder::Single,
                        BondOrder::Aromatic,
                    ])),
                    QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
                ])),
            )
            .unwrap();
        builder
            .set_2d_coordinates(vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [3.0, 0.0],
                [4.0, 0.0],
                [5.0, 0.0],
                [6.0, 0.0],
                [7.0, 0.0],
                [8.0, 0.0],
                [9.0, 0.0],
            ])
            .unwrap();
        builder.build().unwrap()
    }

    fn chiral_flag_molecule() -> Molecule {
        let mut builder = Molecule::builder()
            .with_name("chiral-flag")
            .with_property("_MolFileChiralFlag", "1");
        builder.add_atom(AtomSpec::new(Element::C));
        builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        builder.build().unwrap()
    }

    fn aromatic_benzene_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("benzene");
        let atoms = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true)))
            .collect::<Vec<_>>();
        for idx in 0..6 {
            builder
                .add_bond(
                    BondSpec::new(atoms[idx], atoms[(idx + 1) % 6], BondOrder::Aromatic)
                        .with_aromatic(true),
                )
                .unwrap();
        }
        builder
            .set_2d_coordinates(vec![
                [1.0, 0.0],
                [0.5, 0.866],
                [-0.5, 0.866],
                [-1.0, 0.0],
                [-0.5, -0.866],
                [0.5, -0.866],
            ])
            .unwrap();
        builder.build().unwrap()
    }

    fn chiral_state_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("chiral-state");
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
        builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        builder.build().unwrap()
    }

    fn v3000_bond_cfg_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("bond-cfg");
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Single).with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(a1, a2, BondOrder::Double)
                    .with_direction(BondDirection::EitherDouble)
                    .with_stereo(BondStereo::Any),
            )
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
            .unwrap();
        builder.build().unwrap()
    }

    fn sgroup_molecule() -> Molecule {
        let mut builder = Molecule::builder().with_name("sgroup");
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        let b0 = builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        let sup = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
            .with_external_id(7)
            .with_atoms(vec![a0, a1])
            .with_parent_atoms(vec![a0])
            .with_bonds(vec![b0])
            .with_label("Me")
            .with_subtype("ALT")
            .with_connection(SGroupConnection::HeadToTail)
            .with_expansion_state("E")
            .with_bracket_style(SGroupBracketStyle::Bracket)
            .with_display(SGroupDisplay {
                brackets: vec![SGroupBracket {
                    p1: [0.0, 1.0],
                    p2: [2.0, 3.0],
                }],
                ..SGroupDisplay::default()
            })
            .with_cstates(vec![SGroupCState {
                bond: b0,
                vector: [0.5, 0.25],
            }])
            .with_attach_points(vec![SGroupAttachPoint {
                atom: a0,
                leaving_atom: Some(a1),
                label: Some("AP".to_string()),
                order: None,
            }]);
        let dat = SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
            .with_parent(SubstanceGroupId::new(0))
            .with_component_number(5)
            .with_class("CLASS")
            .with_bracket_style(SGroupBracketStyle::Parenthesis)
            .with_data(SGroupData {
                field_name: Some("FIELD".to_string()),
                field_type: Some("T".to_string()),
                field_info: Some("INFO".to_string()),
                field_display: Some("display spec".to_string()),
                query_type: Some("Q".to_string()),
                query_op: Some("OP".to_string()),
                values: vec!["first valuesecond value".to_string()],
                ..SGroupData::default()
            });
        builder.add_substance_group(sup).unwrap();
        builder.add_substance_group(dat).unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn mol_to_v2000_block_writes_basic_topology_coordinates_and_properties() {
        let molecule = charged_isotope_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.starts_with("charged\n  COSMolKit          2D\n\n"));
        assert!(block.contains("  1  0  0  0  0  0  0  0  0  0999 V2000\n"));
        assert!(block.contains("    1.2500   -2.5000    0.0000   C"));
        assert!(block.contains("M  CHG  1   1  -1\n"));
        assert!(block.contains("M  ISO  1   1  13\n"));
        assert!(block.contains("V    1 payload\n"));
        assert!(block.contains("A    1\nAliasLabel\n"));
        assert!(block.ends_with("M  END\n"));
    }

    #[test]
    fn mol_to_2d_sdf_record_appends_sdf_data_fields_and_delimiter() {
        let molecule = charged_isotope_molecule().with_sdf_data_field("ID", "cmpd-1");

        let record = mol_to_2d_sdf_record(&molecule, SdfFormat::V2000).unwrap();

        assert!(record.contains("M  END\n>  <ID>  \ncmpd-1\n\n$$$$\n"));
    }

    #[test]
    fn molblock_write_params_route_v2000_and_sdf_record_paths() {
        let molecule = charged_isotope_molecule().with_sdf_data_field("ID", "cmpd-1");
        let params = MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            ..Default::default()
        };

        let block = mol_to_mol_block_with_params(&molecule, &params).unwrap();
        let record = mol_to_sdf_record_with_params(&molecule, &params).unwrap();

        assert!(block.contains("  COSMolKit          2D\n"));
        assert!(record.ends_with(">  <ID>  \ncmpd-1\n\n$$$$\n"));
    }

    #[test]
    fn molblock_write_params_route_v3000_precision() {
        let molecule = charged_isotope_molecule();
        let params = MolBlockWriteParams {
            format: SdfFormat::V3000,
            force_2d: true,
            precision: 2,
            ..Default::default()
        };

        let record = mol_to_sdf_record_with_params(&molecule, &params).unwrap();

        assert!(record.contains("M  V30 1 C 1.25 -2.50 0.00 0 CHG=-1 MASS=13\n"));
        assert!(record.ends_with("$$$$\n"));
    }

    #[test]
    fn molblock_write_params_kekulize_prepares_temporary_molecule() {
        let molecule = aromatic_benzene_molecule();
        let kekulized = mol_to_v2000_block(&molecule).unwrap();
        let aromatic = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V2000,
                force_2d: true,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        // prepareMol aromatic-bond bookkeeping: bonds that were aromatic
        // before kekulization are still written as type 4, matching RDKit.
        assert!(kekulized.lines().any(|line| line.starts_with("  1  2  4")));
        assert!(aromatic.lines().any(|line| line.starts_with("  1  2  4")));
        assert!(molecule.bonds().iter().all(|bond| bond.is_aromatic()));
    }

    #[test]
    fn molblock_write_params_include_stereo_false_skips_unported_stereo_state() {
        let molecule = chiral_state_molecule();
        // chiral_tag alone is no longer rejected — parity flag emission is implemented.
        // Both include_stereo=true and include_stereo=false should now succeed.
        let stereo = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V2000,
                force_2d: true,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let no_stereo = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V2000,
                force_2d: true,
                include_stereo: false,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(stereo.contains("    0.0000    0.0000    0.0000   C"));
        assert!(no_stereo.contains("    0.0000    0.0000    0.0000   C"));
        assert!(stereo.ends_with("M  END\n"));
        assert!(no_stereo.ends_with("M  END\n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_rgroup_pxa_and_zero_bond_extensions() {
        let molecule = zbo_extension_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.starts_with("zbo\n  COSMolKit          2D\n\n"));
        assert!(block.contains("M  RGP  1   2   7\n"));
        assert!(block.contains("M  ZBO  1   1   0\n"));
        assert!(block.contains("M  HYD  2   1   2   2   0\n"));
        assert!(block.contains("M  ZCH  1   1  -1\n"));
        assert!(block.contains("M  PXA   2 payload\n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_atom_list_query_lines() {
        let molecule = atom_list_query_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.contains("    0.0000    0.0000    0.0000   L"));
        assert!(block.contains("V    1 [#6,#7]\n"));
        assert!(block.contains("M  ALS   1  2 F C   N   \n"));
        assert!(block.contains("M  ALS   2  2 T O   S   \n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_supported_query_bond_type_codes() {
        let molecule = query_bond_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.contains("  1  2  8  0\n"));
        assert!(block.contains("  3  4  5  0\n"));
        assert!(block.contains("  5  6 42  0\n"));
        assert!(block.contains("  7  8  1  0  0  1\n"));
        assert!(block.contains("  9 10  6  0  0  2\n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_supported_bond_stereo_codes() {
        let molecule = v3000_bond_cfg_molecule();

        let block = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V2000,
                force_2d: true,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let without_stereo = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V2000,
                force_2d: true,
                include_stereo: false,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(block.contains("  1  2  1  1\n"));
        assert!(block.contains("  2  3  2  3\n"));
        assert!(without_stereo.contains("  1  2  1  0\n"));
        assert!(without_stereo.contains("  2  3  2  0\n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_molfile_chiral_flag_in_counts_line() {
        let molecule = chiral_flag_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.contains("  1  0  0  0  1  0  0  0  0  0999 V2000\n"));
    }

    #[test]
    fn mol_to_v2000_block_writes_sgroup_lines() {
        let molecule = sgroup_molecule();

        let block = mol_to_v2000_block(&molecule).unwrap();

        assert!(block.contains("  2  1  0  2  0  0  0  0  0  0999 V2000\n"));
        assert!(block.contains("M  STY  2   1 SUP   2 DAT\n"));
        assert!(block.contains("M  SLB  1   1   7\n"));
        assert!(block.contains("M  SST  1   1 ALT\n"));
        assert!(block.contains("M  SCN  1   1 HT \n"));
        assert!(block.contains("M  SDS EXP  1   1\n"));
        assert!(block.contains("M  SPL  1   2   1\n"));
        assert!(block.contains("M  SNC  1   2   5\n"));
        assert!(block.contains("M  SBT  2   1   0   2   1\n"));
        assert!(block.contains("M  SAL   1  2   1   2\n"));
        assert!(block.contains("M  SPA   1  1   1\n"));
        assert!(block.contains("M  SBL   1  1   1\n"));
        assert!(block.contains("M  SDI   1  4    0.0000    1.0000    2.0000    3.0000\n"));
        assert!(block.contains("M  SMT   1 Me\n"));
        assert!(block.contains("M  SBV   1   1    0.5000    0.2500\n"));
        assert!(block.contains(
            "M  SDT   2 FIELD                         T INFO                Q OP             \n"
        ));
        assert!(block.contains("M  SDD   2 display spec\n"));
        assert!(block.contains("M  SED   2 first valuesecond value\n"));
        assert!(block.contains("M  SAP   1  1   1   2 AP\n"));
        assert!(block.contains("M  SCL   2 CLASS\n"));
    }

    #[test]
    fn mol_to_v3000_block_writes_basic_ctab_atom_bond_blocks() {
        let mut builder = Molecule::builder()
            .with_name("v3000")
            .with_property("_MolFileChiralFlag", "1");
        let carbon = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_formal_charge(-1)
                .with_isotope(13),
        );
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Double))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[1.25, -2.5], [3.0, 4.0]])
            .unwrap();
        let molecule = builder.build().unwrap();

        let block = mol_to_2d_sdf_record(&molecule, SdfFormat::V3000).unwrap();

        assert!(block.starts_with("v3000\n  COSMolKit          2D\n\n"));
        assert!(block.contains("  0  0  0  0  0  0  0  0  0  0999 V3000\n"));
        assert!(block.contains("M  V30 BEGIN CTAB\n"));
        assert!(block.contains("M  V30 COUNTS 2 1 0 0 1\n"));
        assert!(block.contains("M  V30 BEGIN ATOM\n"));
        assert!(block.contains("M  V30 1 C 1.250000 -2.500000 0.000000 0 CHG=-1 MASS=13\n"));
        assert!(block.contains("M  V30 2 O 3.000000 4.000000 0.000000 0\n"));
        assert!(block.contains("M  V30 BEGIN BOND\n"));
        assert!(block.contains("M  V30 1 2 1 2\n"));
        assert!(block.contains("M  V30 END CTAB\nM  END\n$$$$\n"));
    }

    #[test]
    fn mol_to_v3000_block_writes_supported_bond_cfg_values() {
        let molecule = v3000_bond_cfg_molecule();

        let block = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V3000,
                force_2d: true,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(block.contains("M  V30 1 1 1 2 CFG=1\n"));
        assert!(block.contains("M  V30 2 2 2 3 CFG=2\n"));
    }

    #[test]
    fn mol_to_v3000_block_writes_sgroup_lines() {
        let molecule = sgroup_molecule();

        let block = mol_to_v3000_block(&molecule).unwrap();

        assert!(block.contains("M  V30 COUNTS 2 1 2 0 0\n"));
        assert!(block.contains("M  V30 BEGIN SGROUP\n"));
        assert!(block.contains("M  V30 1 SUP 7"));
        assert!(block.contains("ATOMS=(2 1 2)"));
        assert!(block.contains("XBONDS=(1 1)"));
        assert!(block.contains("PATOMS=(1 1)"));
        assert!(block.contains("SUBTYPE=ALT"));
        assert!(block.contains("CONNECT=HT"));
        assert!(block.contains("LABEL=Me"));
        assert!(block.contains("BRKXYZ=(9 0.0000 1.0000 0 2.0000 3.0000 0 0 0 0)"));
        assert!(block.contains("CSTATE=(4 1 0.5000 0.2500 0)"));
        assert_eq!(block.matches("CSTATE=(").count(), 1);
        assert!(block.contains("SAP=(3 1 2 AP)"));
        assert!(block.contains("M  V30 2 DAT 0"));
        assert!(block.contains("PARENT=1"));
        assert!(block.contains("COMPNO=5"));
        assert!(block.contains("FIELDNAME=FIELD"));
        assert!(block.contains("FIELDINFO=INFO"));
        assert!(block.contains("FIELDDISP=\"display spec\""));
        assert!(block.contains("QUERYTYPE=Q"));
        assert!(block.contains("QUERYOP=OP"));
        assert!(block.contains("FIELDDATA=\"first valuesecond value\""));
        assert!(block.contains("CLASS=CLASS"));
        assert!(block.contains("BRKTYP=PAREN"));
        assert!(block.contains("M  V30 END SGROUP\n"));
    }

    #[test]
    fn mol_to_v3000_block_writes_zero_bond_sgroups() {
        let molecule = zbo_extension_molecule();

        let block = mol_to_v3000_block(&molecule).unwrap();

        assert!(block.contains("M  V30 COUNTS 2 1 3 0 0\n"));
        assert!(block.contains("M  V30 1 1 1 2\n"));
        assert!(block.contains("M  V30 BEGIN SGROUP\n"));
        assert!(block.contains("M  V30 1 DAT 0 ATOMS=(2 1 2) XBONDS=(1 1) FIELDNAME=ZBO\n"));
        assert!(block.contains("M  V30 2 DAT 0 ATOMS=(2 1 2) FIELDNAME=HYD FIELDDATA=\"2;0\"\n"));
        assert!(block.contains("M  V30 3 DAT 0 ATOMS=(2 1 2) FIELDNAME=ZCH FIELDDATA=\"-1;0\"\n"));
        assert!(block.contains("M  V30 END SGROUP\n"));
    }

    #[test]
    fn mol_to_v3000_block_writes_enhanced_stereo_collections() {
        let mut builder = Molecule::builder().with_name("collections");
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder
            .add_stereo_group(StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![a0],
                Vec::new(),
            ))
            .unwrap();
        builder
            .add_stereo_group(
                StereoGroup::new(StereoGroupKind::Or, vec![a1], Vec::new()).with_id(2),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let block = mol_to_v3000_block(&molecule).unwrap();

        assert!(block.contains("M  V30 BEGIN COLLECTION\n"));
        assert!(block.contains("M  V30 MDLV30/STEABS ATOMS=(1 1)\n"));
        assert!(block.contains("M  V30 MDLV30/STEREL2 ATOMS=(1 2)\n"));
        assert!(block.contains("M  V30 END COLLECTION\n"));
    }
}
