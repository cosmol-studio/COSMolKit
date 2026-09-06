//! Source-shaped shared state built before Avalon fingerprint feature families.

use crate::FingerprintError;
use crate::chemistry::valence::rdkit_atomic_number_from_symbol;

use super::AvalonFingerprintFlags;
use super::aromaticity::perceive_aromatic_bonds;
use super::daylight_aromaticity::perceive_daylight_aromaticity;
use super::hash::hash_string;
use super::preprocess::{
    Neighbourhood, collect_hydrogen_counts, ring_state, set_ring_size_flags, setup_neighbourhood,
};
use super::reaccs::MoleculeState;
use super::symbols::atom_symbol_match;

const SINGLE: i32 = 1;
const DOUBLE: i32 = 2;
const TRIPLE: i32 = 3;
const AROMATIC: i32 = 4;
pub(super) const USE_DY_AROMATICITY: i32 = 0x0001;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FingerprintPreprocessingState {
    pub neighbours: Vec<Neighbourhood>,
    pub old_bond_types: Vec<i32>,
    pub hydrogen_counts: Vec<i32>,
    pub atom_status: Vec<i32>,
    pub bond_status: Vec<i32>,
    pub degree: Vec<i32>,
    pub carbon_degree: Vec<i32>,
    pub unsaturated: Vec<bool>,
    pub special_neighbours: Vec<i32>,
    pub rare_atom_count: i32,
    pub double_bond_count: i32,
    pub aromatic_bond_count: i32,
    pub fusion_bond_count: i32,
}

fn validate_l_symbol_lists(molecule: &MoleculeState) -> Result<(), FingerprintError> {
    // The Avalon source reads an uninitialized `is_rare` value when an L atom
    // has no matching symbol-list node. COSMolKit fails closed before any
    // query-H or aromaticity mutation while preserving the source behavior for
    // every defined input.
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        if atom.atom_symbol == "L"
            && !molecule
                .symbol_lists
                .iter()
                .any(|symbol_list| symbol_list.atom == atom_index as i32 + 1)
        {
            return Err(FingerprintError::AvalonConversion {
                reason: format!(
                    "Avalon L atom {} has no matching symbol list",
                    atom_index + 1
                ),
            });
        }
    }
    Ok(())
}

pub(super) fn prepare_fingerprint_state(
    molecule: &mut MoleculeState,
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    fpflags: i32,
    exclude_atom: i32,
) -> Result<FingerprintPreprocessingState, FingerprintError> {
    validate_l_symbol_lists(molecule)?;
    // Avalon❗✔️:    nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
    // Avalon❗✔️:    SetupNeighbourhood(mp,nbp,mp->n_atoms);
    let neighbours = setup_neighbourhood(molecule, molecule.atoms.len())?;
    // Avalon❗✔️:    touched_indices = TypeAlloc(mp->n_atoms, int);
    // The traversal-owned touched array is allocated with the traversal unit.
    // Avalon❗✔️:
    // Avalon❗✔️:    old_bond_types = TypeAlloc(mp->n_bonds, int);
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:       old_bond_types[i] = mp->bond_array[i].bond_type;
    let old_bond_types = molecule
        .bonds
        .iter()
        .map(|bond| bond.bond_type)
        .collect::<Vec<_>>();
    // Avalon❗✔️:
    // The complete source H_count allocation and query/non-query branches are
    // reproduced line-by-line in `collect_hydrogen_counts` and its helpers.
    let hydrogen_counts = collect_hydrogen_counts(molecule, &neighbours, as_query)?;
    // Avalon❗✔️:
    // Avalon❗✔️:    /* Do some perception */
    // Avalon❗✔️:    if (fpflags & USE_DY_AROMATICITY)
    // Avalon❗✔️:    {
    if fpflags & USE_DY_AROMATICITY != 0 {
        // Avalon❗✔️:       PerceiveDYAromaticity(mp, nbp);
        perceive_daylight_aromaticity(molecule, &neighbours)?;
        // Avalon❗✔️:    }
        // Avalon❗✔️:    else
        // Avalon❗✔️:    {
    } else {
        // Avalon❗✔️:       PerceiveAromaticBonds(mp);
        perceive_aromatic_bonds(molecule)?;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    degree = TypeAlloc(mp->n_atoms, int);
    let mut degree = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:    cdegree = TypeAlloc(mp->n_atoms, int);
    let mut carbon_degree = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:    nspecial = TypeAlloc(mp->n_atoms, int);
    let mut special_neighbours = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:    unsaturated = TypeAlloc(mp->n_atoms, int);
    let mut unsaturated = vec![false; molecule.atoms.len()];
    // Avalon❗✔️:    bond_status = TypeAlloc(mp->n_bonds, int);
    // Avalon❗✔️:    atom_status = TypeAlloc(mp->n_atoms, int);
    // Avalon❗✔️:    RingState(mp, atom_status, bond_status);
    let (atom_status, bond_status) = ring_state(molecule)?;
    // Avalon❗✔️:    SetRingSizeFlags(mp, 14, nbp);
    set_ring_size_flags(molecule, 14, &neighbours)?;
    // Avalon❗✔️:
    // Avalon❗✔️:    nrare_atoms = 0;
    let mut rare_atom_count = 0_i32;
    // Avalon❗✔️:    /* Set the color property to represent all different atom types */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for atom_index in 0..molecule.atoms.len() {
        // Avalon❗✔️:       unsaturated[i] = FALSE;
        // The calloc-equivalent bool vector is already false.
        let atom_symbol = molecule.atoms[atom_index].atom_symbol.as_str();
        // Avalon❗✔️:       ap->color = StringToInt(periodic_table, ap->atom_symbol);
        let mut color = avalon_atomic_number(atom_symbol);
        // Avalon❗✔️:       if (ap->color <= 1) ap->color = 0;        /* ignore hydrogens */
        if color <= 1 {
            color = 0;
        }
        // Avalon❗✔️:       /* mark special atom types */
        // Avalon❗✔️:       if (ap->color > 115) ap->color = -1;
        if color > 115 {
            color = -1;
        }
        // Avalon❗✔️:       if (0 == strcmp("A", ap->atom_symbol)) ap->color = -1;
        if atom_symbol == "A" {
            color = -1;
        }
        // Avalon❗✔️:       if (0 == strcmp("L", ap->atom_symbol))
        // Avalon❗✔️:       {
        let is_rare = if atom_symbol == "L" {
            // Avalon❗✔️:          ap->color = -1;
            color = -1;
            // Avalon❗✔️:          symbol_lists = mp->symbol_lists;
            let mut list_result = None;
            // Avalon❗✔️:          while (!IsNULL(symbol_lists))
            // Avalon❗✔️:          {
            for symbol_list in &molecule.symbol_lists {
                // Avalon❗✔️:             if (symbol_lists->atom == i+1)
                // Avalon❗✔️:             {
                if symbol_list.atom == atom_index as i32 + 1 {
                    // Avalon❗✔️:                is_rare = TRUE;
                    let mut rare = true;
                    // Avalon❗✔️:                if (!symbol_lists->logic) is_rare = FALSE;
                    if !symbol_list.inclusive {
                        rare = false;
                    }
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("C",symbol_lists->string);
                    rare = rare && !atom_symbol_match("C", &symbol_list.symbols);
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("H",symbol_lists->string);
                    rare = rare && !atom_symbol_match("H", &symbol_list.symbols);
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("O",symbol_lists->string);
                    rare = rare && !atom_symbol_match("O", &symbol_list.symbols);
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("N",symbol_lists->string);
                    rare = rare && !atom_symbol_match("N", &symbol_list.symbols);
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("Cl",symbol_lists->string);
                    rare = rare && !atom_symbol_match("Cl", &symbol_list.symbols);
                    // Avalon❗✔️:                is_rare = is_rare && !AtomSymbolMatch("F",symbol_lists->string);
                    rare = rare && !atom_symbol_match("F", &symbol_list.symbols);
                    // Avalon❗✔️:             }
                    // Preserve the source's last-matching-list-wins order.
                    list_result = Some(rare);
                }
                // Avalon❗✔️:             symbol_lists = symbol_lists->next;
            }
            // Avalon❗✔️:          }
            list_result.ok_or_else(|| FingerprintError::AvalonConversion {
                reason: format!(
                    "Avalon L atom {} has no matching symbol list",
                    atom_index + 1
                ),
            })?
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else
        } else {
            // Avalon❗✔️:          is_rare = ap->color > 0  &&
            // Avalon❗✔️:                    !AtomSymbolMatch(ap->atom_symbol,"C,H,O,N,S,P,Cl,F");
            color > 0 && !atom_symbol_match(atom_symbol, "C,H,O,N,S,P,Cl,F")
        };
        // Avalon❗✔️:       if (is_rare)
        // Avalon❗✔️:          if (exclude_atom != i+1  ||  exclude_atom <= 0) nrare_atoms++;
        if is_rare && (exclude_atom != atom_index as i32 + 1 || exclude_atom <= 0) {
            rare_atom_count += 1;
        }
        // Avalon❗✔️:       // color shortcut atoms if any
        // Avalon❗✔️:       if ((which_bits & USE_SHORTCUT_LABELS)  &&  strcmp(ap->atom_symbol,"R") == 0  &&  ap->atext[0] != '\0')
        // Avalon❗✔️:       {
        if which_bits.contains(AvalonFingerprintFlags::SHORTCUT_LABELS)
            && atom_symbol == "R"
            && !molecule.atoms[atom_index].atext.is_empty()
        {
            // Avalon❗✔️:          ap->color = (0xFFFF00&hash_string(ap->atext)) | 119;
            color = ((hash_string(&molecule.atoms[atom_index].atext) & 0x00ff_ff00) | 119) as i32;
            // Avalon❗✔️:       }
        }
        molecule.atoms[atom_index].color = color;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    ndouble = 0;
    let mut double_bond_count = 0_i32;
    // Avalon❗✔️:    naromatic = 0;
    let mut aromatic_bond_count = 0_i32;
    // Avalon❗✔️:    nfusionb = 0;
    let mut fusion_bond_count = 0_i32;
    // Avalon❗✔️:    /* Set the color property to represent the different bond type classes */
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:    {
    for bond_index in 0..molecule.bonds.len() {
        let first = molecule.bonds[bond_index].atoms[0] as usize - 1;
        let second = molecule.bonds[bond_index].atoms[1] as usize - 1;
        let bond_type = molecule.bonds[bond_index].bond_type;
        // Avalon❗✔️:       if (bp->bond_type == SINGLE)        bp->color = 1;
        let color = if bond_type == SINGLE {
            1
        // Avalon❗✔️:       else if (bp->bond_type == DOUBLE)
        // Avalon❗✔️:       {
        } else if bond_type == DOUBLE {
            // Avalon❗✔️:          bp->color = 2;
            // Avalon❗✔️:          if (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom)
            // Avalon❗✔️:             ndouble++;
            if molecule.bonds[bond_index].atoms[0] != exclude_atom
                && molecule.bonds[bond_index].atoms[1] != exclude_atom
            {
                double_bond_count += 1;
            }
            2
        // Avalon❗✔️:       }
        // Avalon❗✔️:       else if (bp->bond_type == TRIPLE)   bp->color = 3;
        } else if bond_type == TRIPLE {
            3
        // Avalon❗✔️:       else if (bp->bond_type == AROMATIC)
        // Avalon❗✔️:       {
        } else if bond_type == AROMATIC {
            // Avalon❗✔️:          bp->color = 4;
            // Avalon❗✔️:          if (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom)
            // Avalon❗✔️:             naromatic++;
            if molecule.bonds[bond_index].atoms[0] != exclude_atom
                && molecule.bonds[bond_index].atoms[1] != exclude_atom
            {
                aromatic_bond_count += 1;
            }
            4
        // Avalon❗✔️:       }
        // Avalon❗✔️:       else                                bp->color = 0;
        } else {
            0
        };
        molecule.bonds[bond_index].color = color;
        // Avalon❗✔️:       if (bp->color > 1)
        // Avalon❗✔️:       {
        if color > 1 {
            // Avalon❗✔️:          unsaturated[bp->atoms[0]-1] = TRUE;
            unsaturated[first] = true;
            // Avalon❗✔️:          unsaturated[bp->atoms[1]-1] = TRUE;
            unsaturated[second] = true;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:
        // Avalon❗✔️:       /* Count non-hydrogen degree */
        // Avalon❗✔️:       if (mp->atom_array[bp->atoms[0]-1].color != 0  &&
        // Avalon❗✔️:           mp->atom_array[bp->atoms[1]-1].color != 0)
        // Avalon❗✔️:       {
        if molecule.atoms[first].color != 0 && molecule.atoms[second].color != 0 {
            // Avalon❗✔️:          degree[bp->atoms[0]-1]++;
            degree[first] += 1;
            // Avalon❗✔️:          degree[bp->atoms[1]-1]++;
            degree[second] += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       /* Count carbon degree */
        // Avalon❗✔️:       if (bp->bond_type == DOUBLE)
        // Avalon❗✔️:       {
        if bond_type == DOUBLE {
            // Avalon❗✔️:          nspecial[bp->atoms[0]-1]++;
            special_neighbours[first] += 1;
            // Avalon❗✔️:          nspecial[bp->atoms[1]-1]++;
            special_neighbours[second] += 1;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (bp->bond_type == TRIPLE)
            // Avalon❗✔️:       {
        } else if bond_type == TRIPLE {
            // Avalon❗✔️:          nspecial[bp->atoms[0]-1]+=2;
            special_neighbours[first] += 2;
            // Avalon❗✔️:          nspecial[bp->atoms[1]-1]+=2;
            special_neighbours[second] += 2;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (bp->bond_type != SINGLE) continue;
        if bond_type != SINGLE {
            continue;
        }
        // Avalon❗✔️:       if (mp->atom_array[bp->atoms[0]-1].color != 0  &&
        // Avalon❗✔️:           mp->atom_array[bp->atoms[1]-1].color == 6)
        // Avalon❗✔️:          cdegree[bp->atoms[0]-1]++;
        if molecule.atoms[first].color != 0 && molecule.atoms[second].color == 6 {
            carbon_degree[first] += 1;
        }
        // Avalon❗✔️:       if (mp->atom_array[bp->atoms[1]-1].color != 0  &&
        // Avalon❗✔️:           mp->atom_array[bp->atoms[0]-1].color == 6)
        // Avalon❗✔️:          cdegree[bp->atoms[1]-1]++;
        if molecule.atoms[second].color != 0 && molecule.atoms[first].color == 6 {
            carbon_degree[second] += 1;
        }
        // Avalon❗✔️:
        // Avalon❗✔️:       if (exclude_atom <= 0  ||
        // Avalon❗✔️:           (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom))
        // Avalon❗✔️:          if (atom_status[bp->atoms[0]-1] > 2  &&
        // Avalon❗✔️:              atom_status[bp->atoms[1]-1] > 2) // ring fusion
        // Avalon❗✔️:             nfusionb++;
        if (exclude_atom <= 0
            || (molecule.bonds[bond_index].atoms[0] != exclude_atom
                && molecule.bonds[bond_index].atoms[1] != exclude_atom))
            && atom_status[first] > 2
            && atom_status[second] > 2
        {
            fusion_bond_count += 1;
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    /* ignore special atom types for further processing */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       if (ap->color < 0) ap->color = 0;
    for atom in &mut molecule.atoms {
        if atom.color < 0 {
            atom.color = 0;
        }
    }

    Ok(FingerprintPreprocessingState {
        neighbours,
        old_bond_types,
        hydrogen_counts,
        atom_status,
        bond_status,
        degree,
        carbon_degree,
        unsaturated,
        special_neighbours,
        rare_atom_count,
        double_bond_count,
        aromatic_bond_count,
        fusion_bond_count,
    })
}

pub(super) fn restore_fingerprint_bond_types(
    molecule: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
) -> Result<(), FingerprintError> {
    if molecule.bonds.len() != state.old_bond_types.len() {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon fingerprint traversal changed the bond-table length".to_string(),
        });
    }
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:       mp->bond_array[i].bond_type = old_bond_types[i];
    for (bond, &old_bond_type) in molecule.bonds.iter_mut().zip(&state.old_bond_types) {
        bond.bond_type = old_bond_type;
    }
    // Avalon❗✔️:    MyFree((char *)old_bond_types);
    // Rust releases the owned snapshot with `state`.
    Ok(())
}

pub(super) fn with_prepared_fingerprint_state<T>(
    molecule: &mut MoleculeState,
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    fpflags: i32,
    exclude_atom: i32,
    operation: impl FnOnce(
        &mut MoleculeState,
        &FingerprintPreprocessingState,
    ) -> Result<T, FingerprintError>,
) -> Result<T, FingerprintError> {
    let state = prepare_fingerprint_state(molecule, which_bits, as_query, fpflags, exclude_atom)?;
    let result = operation(molecule, &state);
    restore_fingerprint_bond_types(molecule, &state)?;
    result
}

pub(super) fn avalon_atomic_number(symbol: &str) -> i32 {
    // Avalon❗✔️: string_int_table periodic_table[] =
    // Avalon❗✔️:    {
    // Avalon❗✔️:       {"H",             1},
    // Avalon❗✔️:       {"D",             1},
    // Avalon❗✔️:       {"T",             1},
    // Avalon❗✔️:       {"He",            2},
    // Avalon❗✔️:       {"Li",            3},
    // Avalon❗✔️:       {"Be",            4},
    // Avalon❗✔️:       {"B",             5},
    // Avalon❗✔️:       {"C",             6},
    // Avalon❗✔️:       {"N",             7},
    // Avalon❗✔️:       {"O",             8},
    // Avalon❗✔️:       {"F",             9},
    // Avalon❗✔️:       {"Ne",           10},
    // Avalon❗✔️:       {"Na",           11},
    // Avalon❗✔️:       {"Mg",           12},
    // Avalon❗✔️:       {"Al",           13},
    // Avalon❗✔️:       {"Si",           14},
    // Avalon❗✔️:       {"P",            15},
    // Avalon❗✔️:       {"S",            16},
    // Avalon❗✔️:       {"Cl",           17},
    // Avalon❗✔️:       {"Ar",           18},
    // Avalon❗✔️:       {"K",            19},
    // Avalon❗✔️:       {"Ca",           20},
    // Avalon❗✔️:       {"Sc",           21},
    // Avalon❗✔️:       {"Ti",           22},
    // Avalon❗✔️:       {"V",            23},
    // Avalon❗✔️:       {"Cr",           24},
    // Avalon❗✔️:       {"Mn",           25},
    // Avalon❗✔️:       {"Fe",           26},
    // Avalon❗✔️:       {"Co",           27},
    // Avalon❗✔️:       {"Ni",           28},
    // Avalon❗✔️:       {"Cu",           29},
    // Avalon❗✔️:       {"Zn",           30},
    // Avalon❗✔️:       {"Ga",           31},
    // Avalon❗✔️:       {"Ge",           32},
    // Avalon❗✔️:       {"As",           33},
    // Avalon❗✔️:       {"Se",           34},
    // Avalon❗✔️:       {"Br",           35},
    // Avalon❗✔️:       {"Kr",           36},
    // Avalon❗✔️:       {"Rb",           37},
    // Avalon❗✔️:       {"Sr",           38},
    // Avalon❗✔️:       {"Y",            39},
    // Avalon❗✔️:       {"Zr",           40},
    // Avalon❗✔️:       {"Nb",           41},
    // Avalon❗✔️:       {"Mo",           42},
    // Avalon❗✔️:       {"Tc",           43},
    // Avalon❗✔️:       {"Ru",           44},
    // Avalon❗✔️:       {"Rh",           45},
    // Avalon❗✔️:       {"Pd",           46},
    // Avalon❗✔️:       {"Ag",           47},
    // Avalon❗✔️:       {"Cd",           48},
    // Avalon❗✔️:       {"In",           49},
    // Avalon❗✔️:       {"Sn",           50},
    // Avalon❗✔️:       {"Sb",           51},
    // Avalon❗✔️:       {"Te",           52},
    // Avalon❗✔️:       {"I",            53},
    // Avalon❗✔️:       {"Xe",           54},
    // Avalon❗✔️:       {"Cs",           55},
    // Avalon❗✔️:       {"Ba",           56},
    // Avalon❗✔️:       {"La",           57},
    // Avalon❗✔️:       {"Ce",           58},
    // Avalon❗✔️:       {"Pr",           59},
    // Avalon❗✔️:       {"Nd",           60},
    // Avalon❗✔️:       {"Pm",           61},
    // Avalon❗✔️:       {"Sm",           62},
    // Avalon❗✔️:       {"Eu",           63},
    // Avalon❗✔️:       {"Gd",           64},
    // Avalon❗✔️:       {"Tb",           65},
    // Avalon❗✔️:       {"Dy",           66},
    // Avalon❗✔️:       {"Ho",           67},
    // Avalon❗✔️:       {"Er",           68},
    // Avalon❗✔️:       {"Tm",           69},
    // Avalon❗✔️:       {"Yb",           70},
    // Avalon❗✔️:       {"Lu",           71},
    // Avalon❗✔️:       {"Hf",           72},
    // Avalon❗✔️:       {"Ta",           73},
    // Avalon❗✔️:       {"W",            74},
    // Avalon❗✔️:       {"Re",           75},
    // Avalon❗✔️:       {"Os",           76},
    // Avalon❗✔️:       {"Ir",           77},
    // Avalon❗✔️:       {"Pt",           78},
    // Avalon❗✔️:       {"Au",           79},
    // Avalon❗✔️:       {"Hg",           80},
    // Avalon❗✔️:       {"Tl",           81},
    // Avalon❗✔️:       {"Pb",           82},
    // Avalon❗✔️:       {"Bi",           83},
    // Avalon❗✔️:       {"Po",           84},
    // Avalon❗✔️:       {"At",           85},
    // Avalon❗✔️:       {"Rn",           86},
    // Avalon❗✔️:       {"Fr",           87},
    // Avalon❗✔️:       {"Ra",           88},
    // Avalon❗✔️:       {"Ac",           89},
    // Avalon❗✔️:       {"Th",           90},
    // Avalon❗✔️:       {"Pa",           91},
    // Avalon❗✔️:       {"U",            92},
    // Avalon❗✔️:       {"Np",           93},
    // Avalon❗✔️:       {"Pu",           94},
    // Avalon❗✔️:       {"X",           120},
    // Avalon❗✔️:       {"Q",           121},
    // Avalon❗✔️:       {"M",           122},
    // Avalon❗✔️:       {"R",           123},
    // Avalon❗✔️:       {"A",           124},
    // Avalon❗✔️:       {(char *)NULL,    0}
    // Avalon❗✔️:    };
    match symbol {
        "D" | "T" => 1,
        "X" => 120,
        "Q" => 121,
        "M" => 122,
        "R" => 123,
        "A" => 124,
        _ => rdkit_atomic_number_from_symbol(symbol)
            .filter(|&number| number <= 94)
            .map_or(0, i32::from),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond, SymbolList};

    fn atom(symbol: &str) -> Atom {
        Atom {
            atom_symbol: symbol.to_string(),
            ..Atom::default()
        }
    }

    fn bond(atoms: [i32; 2], bond_type: i32) -> Bond {
        Bond {
            atoms,
            bond_type,
            ..Bond::default()
        }
    }

    fn colors(molecule: &MoleculeState) -> (Vec<i32>, Vec<i32>) {
        (
            molecule.atoms.iter().map(|atom| atom.color).collect(),
            molecule.bonds.iter().map(|bond| bond.color).collect(),
        )
    }

    #[test]
    fn mixed_atom_bond_state_matches_native_intermediate_arrays() {
        let mut molecule = MoleculeState {
            atoms: ["C", "N", "O", "H", "Fe"].map(atom).to_vec(),
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([1, 3], DOUBLE),
                bond([1, 4], SINGLE),
                bond([2, 5], TRIPLE),
            ],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();

        assert_eq!(colors(&molecule), (vec![6, 7, 8, 0, 26], vec![1, 2, 1, 3]));
        assert_eq!(state.degree, vec![2, 2, 1, 0, 1]);
        assert_eq!(state.carbon_degree, vec![0, 1, 0, 0, 0]);
        assert_eq!(state.special_neighbours, vec![1, 2, 1, 0, 2]);
        assert_eq!(state.unsaturated, vec![true, true, true, false, true]);
        assert_eq!(state.atom_status, vec![0; 5]);
        assert_eq!(state.bond_status, vec![0; 4]);
        assert!(molecule.atoms.iter().all(|atom| atom.rsize_flags == 0));
        assert!(molecule.bonds.iter().all(|bond| bond.rsize_flags == 0));
        assert_eq!(state.rare_atom_count, 1);
        assert_eq!(state.double_bond_count, 1);
        assert_eq!(state.aromatic_bond_count, 0);
        assert_eq!(state.fusion_bond_count, 0);
    }

    #[test]
    fn exclude_atom_changes_only_source_counters() {
        let mut molecule = MoleculeState {
            atoms: ["C", "N", "O", "H", "Fe"].map(atom).to_vec(),
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([1, 3], DOUBLE),
                bond([1, 4], SINGLE),
                bond([2, 5], TRIPLE),
            ],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            USE_DY_AROMATICITY,
            5,
        )
        .unwrap();

        assert_eq!(state.rare_atom_count, 0);
        assert_eq!(state.degree, vec![2, 2, 1, 0, 1]);
        assert_eq!(state.special_neighbours, vec![1, 2, 1, 0, 2]);
    }

    #[test]
    fn fused_ring_state_matches_native_counters() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 4],
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([2, 3], DOUBLE),
                bond([3, 4], SINGLE),
                bond([4, 1], DOUBLE),
                bond([1, 3], SINGLE),
            ],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();

        assert_eq!(colors(&molecule), (vec![6; 4], vec![1, 2, 1, 2, 1]));
        assert_eq!(state.degree, vec![3, 2, 3, 2]);
        assert_eq!(state.carbon_degree, vec![2, 1, 2, 1]);
        assert_eq!(state.special_neighbours, vec![1; 4]);
        assert_eq!(state.unsaturated, vec![true; 4]);
        assert_eq!(state.atom_status, vec![3, 2, 3, 2]);
        assert_eq!(state.bond_status, vec![1, 1, 1, 1, 2]);
        assert_eq!(
            molecule
                .atoms
                .iter()
                .map(|atom| atom.rsize_flags)
                .collect::<Vec<_>>(),
            vec![25; 4]
        );
        assert_eq!(
            molecule
                .bonds
                .iter()
                .map(|bond| bond.rsize_flags)
                .collect::<Vec<_>>(),
            vec![25, 25, 25, 25, 9]
        );
        assert_eq!(state.double_bond_count, 2);
        assert_eq!(state.fusion_bond_count, 1);
    }

    #[test]
    fn symbol_lists_and_shortcut_color_match_native_source_order() {
        let mut shortcut = atom("R");
        shortcut.atext = "ALA".to_string();
        let mut molecule = MoleculeState {
            atoms: vec![atom("L"), shortcut, atom("L")],
            symbol_lists: vec![
                SymbolList {
                    atom: 1,
                    inclusive: true,
                    symbols: "Br,I".to_string(),
                },
                SymbolList {
                    atom: 3,
                    inclusive: true,
                    symbols: "C,N".to_string(),
                },
            ],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::SHORTCUT_LABELS,
            false,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();

        assert_eq!(colors(&molecule).0, vec![0, 14_602_103, 0]);
        assert_eq!(state.rare_atom_count, 1);
    }

    #[test]
    fn missing_l_symbol_list_is_structured_before_source_uninitialized_read() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("L")],
            ..MoleculeState::default()
        };

        assert!(matches!(
            prepare_fingerprint_state(
                &mut molecule,
                AvalonFingerprintFlags::from_bits_retain(0),
                true,
                USE_DY_AROMATICITY,
                0,
            ),
            Err(FingerprintError::AvalonConversion { reason })
                if reason == "Avalon L atom 1 has no matching symbol list"
        ));
    }

    #[test]
    fn missing_l_symbol_list_fails_before_aromaticity_changes_bonds() {
        let original = vec![SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut molecule = MoleculeState {
            atoms: ["L", "C", "C", "C", "C", "C"].map(atom).to_vec(),
            bonds: [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
                .into_iter()
                .zip(original.iter().copied())
                .map(|(atoms, bond_type)| bond(atoms, bond_type))
                .collect(),
            ..MoleculeState::default()
        };

        assert!(matches!(
            prepare_fingerprint_state(
                &mut molecule,
                AvalonFingerprintFlags::from_bits_retain(0),
                false,
                0,
                0,
            ),
            Err(FingerprintError::AvalonConversion { reason })
                if reason == "Avalon L atom 1 has no matching symbol list"
        ));
        assert_eq!(
            molecule
                .bonds
                .iter()
                .map(|bond| bond.bond_type)
                .collect::<Vec<_>>(),
            original
        );
    }

    #[test]
    fn query_mode_uses_query_h_count_before_shared_color_state() {
        let mut nitrogen = atom("N");
        nitrogen.query_h_count = 4;
        let mut molecule = MoleculeState {
            atoms: vec![nitrogen],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            true,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();

        assert_eq!(state.hydrogen_counts, vec![0, 3]);
        assert_eq!(colors(&molecule).0, vec![7]);
    }

    #[test]
    fn standard_aromaticity_branch_keeps_the_source_bond_type_snapshot() {
        let original = vec![SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 6],
            bonds: [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
                .into_iter()
                .zip(original.iter().copied())
                .map(|(atoms, bond_type)| bond(atoms, bond_type))
                .collect(),
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            0,
            0,
        )
        .unwrap();

        assert_eq!(state.old_bond_types, original);
        assert!(
            molecule
                .bonds
                .iter()
                .all(|bond| bond.bond_type == AROMATIC && bond.color == 4)
        );
        assert_eq!(state.aromatic_bond_count, 6);
    }

    #[test]
    fn avalon_periodic_table_stops_at_pu_and_keeps_source_pseudo_values() {
        assert_eq!(avalon_atomic_number("H"), 1);
        assert_eq!(avalon_atomic_number("D"), 1);
        assert_eq!(avalon_atomic_number("T"), 1);
        assert_eq!(avalon_atomic_number("He"), 2);
        assert_eq!(avalon_atomic_number("Pu"), 94);
        assert_eq!(avalon_atomic_number("Am"), 0);
        assert_eq!(avalon_atomic_number("X"), 120);
        assert_eq!(avalon_atomic_number("Q"), 121);
        assert_eq!(avalon_atomic_number("M"), 122);
        assert_eq!(avalon_atomic_number("R"), 123);
        assert_eq!(avalon_atomic_number("A"), 124);
        assert_eq!(avalon_atomic_number("L"), 0);
    }

    #[test]
    fn last_matching_symbol_list_controls_rare_atom_count() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("L")],
            symbol_lists: vec![
                SymbolList {
                    atom: 1,
                    inclusive: true,
                    symbols: "Br,I".to_string(),
                },
                SymbolList {
                    atom: 1,
                    inclusive: true,
                    symbols: "C,N".to_string(),
                },
            ],
            ..MoleculeState::default()
        };

        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            true,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();
        assert_eq!(state.rare_atom_count, 0);

        molecule.symbol_lists.reverse();
        let state = prepare_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            true,
            USE_DY_AROMATICITY,
            0,
        )
        .unwrap();
        assert_eq!(state.rare_atom_count, 1);
    }

    #[test]
    fn complete_preprocessing_lifecycle_restores_bond_types_after_error() {
        let original = vec![SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 6],
            bonds: [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
                .into_iter()
                .zip(original.iter().copied())
                .map(|(atoms, bond_type)| bond(atoms, bond_type))
                .collect(),
            ..MoleculeState::default()
        };

        let result = with_prepared_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            0,
            0,
            |working, _state| {
                assert!(working.bonds.iter().all(|bond| bond.bond_type == AROMATIC));
                Err::<(), _>(FingerprintError::InvalidArguments {
                    reason: "forced traversal error",
                })
            },
        );

        assert!(matches!(
            result,
            Err(FingerprintError::InvalidArguments { .. })
        ));
        assert_eq!(
            molecule
                .bonds
                .iter()
                .map(|bond| bond.bond_type)
                .collect::<Vec<_>>(),
            original
        );
    }

    #[test]
    fn complete_preprocessing_lifecycle_restores_bond_types_after_success() {
        let original = vec![SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 6],
            bonds: [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
                .into_iter()
                .zip(original.iter().copied())
                .map(|(atoms, bond_type)| bond(atoms, bond_type))
                .collect(),
            ..MoleculeState::default()
        };

        let aromatic_count = with_prepared_fingerprint_state(
            &mut molecule,
            AvalonFingerprintFlags::from_bits_retain(0),
            false,
            0,
            0,
            |working, state| {
                assert!(working.bonds.iter().all(|bond| bond.bond_type == AROMATIC));
                Ok::<_, FingerprintError>(state.aromatic_bond_count)
            },
        )
        .unwrap();

        assert_eq!(aromatic_count, 6);
        assert_eq!(
            molecule
                .bonds
                .iter()
                .map(|bond| bond.bond_type)
                .collect::<Vec<_>>(),
            original
        );
    }
}
