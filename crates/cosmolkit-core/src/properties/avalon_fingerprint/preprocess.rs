//! Source-shaped preprocessing state shared by Avalon fingerprint families.

use crate::FingerprintError;

use super::reaccs::MoleculeState;
use super::rings::{combine_rings, ring_list};

const MAX_NEIGHBOURS: usize = 20;
const NONE: i32 = 0;
const DOUBLET: i32 = 2;
const ZERO_COUNT: i32 = 1;
const SINGLE: i32 = 1;
const DOUBLE: i32 = 2;
const TRIPLE: i32 = 3;
const AROMATIC: i32 = 4;
const SUB_AS_IS: i32 = -2;

#[derive(Debug, Clone, Copy)]
struct ValenceEntry {
    atom_type: &'static str,
    from_valence: i32,
    to_valence: i32,
    step_valence: i32,
    lone_pairs: i32,
    pair_deficit: i32,
}

// Avalon❗✔️: struct valence_table_entry      /* This is the table of possible valences */
// Avalon❗✔️:    {                            /* of the chemical elements. Elements not */
// Avalon❗✔️:       char *atom_type;          /* in the table are assumed to have no    */
// Avalon❗✔️:       int   from_valence,       /* valence constraints at all.            */
// Avalon❗✔️:             to_valence,
// Avalon❗✔️:             step_valence,
// Avalon❗✔️:             lone_pairs,
// Avalon❗✔️: 	    pair_deficit;
// Avalon❗✔️:    } valence_table[] =  {
// Avalon❗✔️:                         {"C",  4, 4, 1},        // make common elements hit first
// Avalon❗✔️:                         {"H",  -1, 1, 2},
// Avalon❗✔️:                         {"N",  3, 5, 2, 1},
// Avalon❗✔️:                         {"O",  2, 2, 1, 2},
// Avalon❗✔️:                         {"Cl", 1, 7, 2},
// Avalon❗✔️:                         {"P",  3, 5, 2, 1},
// Avalon❗✔️:                         {"S",  2, 6, 2, 2},
// Avalon❗✔️:                         {"F",  1, 1, 1, 1},
// Avalon❗✔️:
// Avalon❗✔️:                         {"H",  -1, 1, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         {"Li", -1, 1, 2},
// Avalon❗✔️:                         {"Na", -1, 1, 2},
// Avalon❗✔️:                         {"K",  -1, 1, 2},
// Avalon❗✔️:                         {"Rb", -1, 1, 2},
// Avalon❗✔️:                         {"Cs", -1, 1, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         {"Be", -2, 2, 2},
// Avalon❗✔️:                         {"Mg", -2, 2, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         {"B",   3, 3, 1, 0, 1},
// Avalon❗✔️:
// Avalon❗✔️:                         {"Al", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Ga", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"In", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Tl", -3, 3, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         // {"C",  4, 4, 1},
// Avalon❗✔️:                         {"Si", 4, 4, 1},
// Avalon❗✔️:
// Avalon❗✔️:                         // {"N",  3, 5, 2, 1},
// Avalon❗✔️:                         // {"P",  3, 5, 2, 1},
// Avalon❗✔️:                         {"As", 3, 5, 2},
// Avalon❗✔️:                         {"Sb", 3, 5, 2},
// Avalon❗✔️:                         {"Bi", 3, 5, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         // {"O",  2, 2, 1, 2},
// Avalon❗✔️:                         // {"S",  2, 6, 2, 2},
// Avalon❗✔️:                         {"Se", 2, 6, 2},
// Avalon❗✔️:                         {"Te", 2, 6, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         {"La", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Ce", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Pr", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Nd", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Pm", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Sm", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Eu", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Gd", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Tb", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Dy", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Ho", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Er", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Tm", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Yb", -3, 3, 2, 0, 1},
// Avalon❗✔️:                         {"Lu", -3, 3, 2, 0, 1},
// Avalon❗✔️:
// Avalon❗✔️:                         // {"F",  1, 1, 1, 1},
// Avalon❗✔️:                         // {"Cl", 1, 7, 2},
// Avalon❗✔️:                         {"Br", 1, 7, 2},
// Avalon❗✔️:                         {"I",  1, 7, 2},
// Avalon❗✔️:
// Avalon❗✔️:                         {(char *)NULL, 0}};
const VALENCE_TABLE: &[ValenceEntry] = &[
    valence("C", 4, 4, 1, 0, 0),
    valence("H", -1, 1, 2, 0, 0),
    valence("N", 3, 5, 2, 1, 0),
    valence("O", 2, 2, 1, 2, 0),
    valence("Cl", 1, 7, 2, 0, 0),
    valence("P", 3, 5, 2, 1, 0),
    valence("S", 2, 6, 2, 2, 0),
    valence("F", 1, 1, 1, 1, 0),
    valence("H", -1, 1, 2, 0, 0),
    valence("Li", -1, 1, 2, 0, 0),
    valence("Na", -1, 1, 2, 0, 0),
    valence("K", -1, 1, 2, 0, 0),
    valence("Rb", -1, 1, 2, 0, 0),
    valence("Cs", -1, 1, 2, 0, 0),
    valence("Be", -2, 2, 2, 0, 0),
    valence("Mg", -2, 2, 2, 0, 0),
    valence("B", 3, 3, 1, 0, 1),
    valence("Al", -3, 3, 2, 0, 1),
    valence("Ga", -3, 3, 2, 0, 1),
    valence("In", -3, 3, 2, 0, 1),
    valence("Tl", -3, 3, 2, 0, 0),
    valence("Si", 4, 4, 1, 0, 0),
    valence("As", 3, 5, 2, 0, 0),
    valence("Sb", 3, 5, 2, 0, 0),
    valence("Bi", 3, 5, 2, 0, 0),
    valence("Se", 2, 6, 2, 0, 0),
    valence("Te", 2, 6, 2, 0, 0),
    valence("La", -3, 3, 2, 0, 1),
    valence("Ce", -3, 3, 2, 0, 1),
    valence("Pr", -3, 3, 2, 0, 1),
    valence("Nd", -3, 3, 2, 0, 1),
    valence("Pm", -3, 3, 2, 0, 1),
    valence("Sm", -3, 3, 2, 0, 1),
    valence("Eu", -3, 3, 2, 0, 1),
    valence("Gd", -3, 3, 2, 0, 1),
    valence("Tb", -3, 3, 2, 0, 1),
    valence("Dy", -3, 3, 2, 0, 1),
    valence("Ho", -3, 3, 2, 0, 1),
    valence("Er", -3, 3, 2, 0, 1),
    valence("Tm", -3, 3, 2, 0, 1),
    valence("Yb", -3, 3, 2, 0, 1),
    valence("Lu", -3, 3, 2, 0, 1),
    valence("Br", 1, 7, 2, 0, 0),
    valence("I", 1, 7, 2, 0, 0),
];

const fn valence(
    atom_type: &'static str,
    from_valence: i32,
    to_valence: i32,
    step_valence: i32,
    lone_pairs: i32,
    pair_deficit: i32,
) -> ValenceEntry {
    ValenceEntry {
        atom_type,
        from_valence,
        to_valence,
        step_valence,
        lone_pairs,
        pair_deficit,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Neighbourhood {
    n_ligands: usize,
    atoms: [usize; MAX_NEIGHBOURS],
    bonds: [usize; MAX_NEIGHBOURS],
}

impl Default for Neighbourhood {
    fn default() -> Self {
        Self {
            n_ligands: 0,
            atoms: [0; MAX_NEIGHBOURS],
            bonds: [0; MAX_NEIGHBOURS],
        }
    }
}

impl Neighbourhood {
    pub(super) fn atoms(&self) -> &[usize] {
        &self.atoms[..self.n_ligands]
    }

    pub(super) fn bonds(&self) -> &[usize] {
        &self.bonds[..self.n_ligands]
    }

    fn push(&mut self, atom: usize, bond: usize, owner: usize) -> Result<(), FingerprintError> {
        if self.n_ligands == MAX_NEIGHBOURS {
            return Err(FingerprintError::AvalonConversion {
                reason: format!("Avalon neighbourhood capacity exceeded at atom {}", owner + 1),
            });
        }
        self.atoms[self.n_ligands] = atom;
        self.bonds[self.n_ligands] = bond;
        self.n_ligands += 1;
        Ok(())
    }
}

pub(super) fn setup_neighbourhood(
    molecule: &MoleculeState,
    nlimit: usize,
) -> Result<Vec<Neighbourhood>, FingerprintError> {
    // Avalon❗✔️: int SetupNeighbourhood(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                        neighbourhood_t *neighbour_array,
    // Avalon❗✔️:                        int nlimit)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Computes the array of atom and bond neighbourhoods of the
    // Avalon❗✔️:  * atoms in *mp which have an atom number less than nlimit.
    // Avalon❗✔️:  *
    // Avalon❗✔️:  * Returns FALSE in case of failure.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    int at0, at1;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)        /* setup neighbourhood */
    // Avalon❗✔️:       neighbour_array[i].n_ligands = 0;
    let mut neighbours = vec![Neighbourhood::default(); molecule.atoms.len()];
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:    {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Avalon❗✔️:       at0 = mp->bond_array[i].atoms[0]-1;
        // Avalon❗✔️:       at1 = mp->bond_array[i].atoms[1]-1;
        let at0 = source_atom_index(bond.atoms[0], molecule.atoms.len(), bond_index)?;
        let at1 = source_atom_index(bond.atoms[1], molecule.atoms.len(), bond_index)?;
        // Avalon❗✔️:       if (at0 >= nlimit  ||  at1 >= nlimit) continue;
        if at0 >= nlimit || at1 >= nlimit {
            continue;
        }
        // Avalon❗✔️:       neighbour_array[at0].atoms[neighbour_array[at0].n_ligands] = at1;
        // Avalon❗✔️:       neighbour_array[at0].bonds[neighbour_array[at0].n_ligands] = i;
        // Avalon❗✔️:       neighbour_array[at0].n_ligands++;
        // Avalon❗✔️:       if (neighbour_array[at0].n_ligands > MAXNEIGHBOURS)
        // Avalon❗✔️:       {
        // Avalon❗✔️:          fprintf(stderr,"Too many ligands at atom %d\n", at0);
        // Avalon❗✔️:          ShowMessageI("Too many neighbours at atom %d\n",
        // Avalon❗✔️:                       "SetupNeighbourhood",
        // Avalon❗✔️:                       at0+1);
        // Avalon❗✔️:          return FALSE;
        // Avalon❗✔️:       }
        // The source checks only after writing element 21 of a 20-element
        // array. Rust rejects that unmodelled undefined-behavior boundary
        // before the write while preserving all valid source states.
        neighbours[at0].push(at1, bond_index, at0)?;
        // Avalon❗✔️:       neighbour_array[at1].atoms[neighbour_array[at1].n_ligands] = at0;
        // Avalon❗✔️:       neighbour_array[at1].bonds[neighbour_array[at1].n_ligands] = i;
        // Avalon❗✔️:       neighbour_array[at1].n_ligands++;
        // Avalon❗✔️:       if (neighbour_array[at1].n_ligands > MAXNEIGHBOURS)
        // Avalon❗✔️:       {
        // Avalon❗✔️:          fprintf(stderr,"Too many ligands at atom %d\n", at1);
        // Avalon❗✔️:          ShowMessageI("Too many neighbours at atom %d\n",
        // Avalon❗✔️:                       "SetupNeighbourhood",
        // Avalon❗✔️:                       at1+1);
        // Avalon❗✔️:          return FALSE;
        // Avalon❗✔️:       }
        neighbours[at1].push(at0, bond_index, at1)?;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    return TRUE;
    // Avalon❗✔️: }
    Ok(neighbours)
}

fn source_atom_index(source_index: i32, atom_count: usize, bond_index: usize) -> Result<usize, FingerprintError> {
    let index = source_index
        .checked_sub(1)
        .and_then(|value| usize::try_from(value).ok());
    index
        .filter(|&value| value < atom_count)
        .ok_or_else(|| FingerprintError::AvalonConversion {
            reason: format!(
                "Avalon bond {} references invalid atom {}",
                bond_index + 1,
                source_index
            ),
        })
}

fn implicit_hydrogens(
    symbol: &str,
    nsingle: i32,
    naromatic: i32,
    ndouble: i32,
    ntriple: i32,
    radical: i32,
    charge: i32,
) -> i32 {
    // Avalon❗✔️: int
    // Avalon❗✔️: ImplicitHydrogens(char *symbol,
    // Avalon❗✔️:                   int   nsingle,
    // Avalon❗✔️:                   int   naromatic,
    // Avalon❗✔️:                   int   ndouble,
    // Avalon❗✔️:                   int   ntriple,
    // Avalon❗✔️:                   int   radical,
    // Avalon❗✔️:                   int   charge)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Computes the number of implicit hydrogens attached to an atom of type
    // Avalon❗✔️:  * symbol with nsingle single bonds, naromatic aromatic bonds, ndouble
    // Avalon❗✔️:  * double bonds, ntriple triple bonds, and radical and charge state
    // Avalon❗✔️:  * radical and charge, resp.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i,val,h;
    // Avalon❗✔️:    int bond_electrons;
    // Avalon❗✔️:
    // Avalon❗✔️:    bond_electrons = nsingle+2*ndouble+3*ntriple;
    let mut bond_electrons = nsingle + 2 * ndouble + 3 * ntriple;
    // Avalon❗✔️:    if (radical) bond_electrons++;
    if radical != 0 {
        bond_electrons += 1;
    }
    // Avalon❗✔️:    switch (naromatic)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       case 0: break;
    // Avalon❗✔️:       case 1: /* one aromatic bond could be an error */
    // Avalon❗✔️:               bond_electrons+=2;
    // Avalon❗✔️:               break;
    // Avalon❗✔️:       case 2: bond_electrons+=3;
    // Avalon❗✔️:               break;
    // Avalon❗✔️:       case 3: bond_electrons+=4;
    // Avalon❗✔️:               break;
    // Avalon❗✔️:       default:ShowMessageI("atom with %d aromatic bonds",
    // Avalon❗✔️:                            "ImplicitHydrogens",
    // Avalon❗✔️:                            naromatic);
    // Avalon❗✔️:               bond_electrons += naromatic+1;
    // Avalon❗✔️:               break;
    // Avalon❗✔️:    }
    bond_electrons += match naromatic {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 4,
        count => count + 1,
    };

    // Avalon❗✔️:    for (i=0; valence_table[i].atom_type!=(char *)NULL; i++)
    // Avalon❗✔️:       if (0 == strcmp(valence_table[i].atom_type,symbol))
    // Avalon❗✔️:       {
    for entry in VALENCE_TABLE.iter().filter(|entry| entry.atom_type == symbol) {
        // Avalon❗✔️:          if (charge == 0)       /* Easy case */
        // Avalon❗✔️:          {
        if charge == 0 {
            // Avalon❗✔️:             for (val=valence_table[i].from_valence;
            // Avalon❗✔️:                  val<=valence_table[i].to_valence;
            // Avalon❗✔️:                  val+=valence_table[i].step_valence)
            // Avalon❗✔️:                if (0 <= (h=val-bond_electrons)) return(h);
            let mut value = entry.from_valence;
            while value <= entry.to_valence {
                let hydrogens = value - bond_electrons;
                if hydrogens >= 0 {
                    return hydrogens;
                }
                value += entry.step_valence;
            }
            // Avalon❗✔️:          }
            // Avalon❗✔️:          else if (charge > 0)
            // Avalon❗✔️:          {
        } else if charge > 0 {
            // Avalon❗✔️:             for (val=valence_table[i].from_valence;
            // Avalon❗✔️:                  val<=valence_table[i].to_valence;
            // Avalon❗✔️:                  val+=valence_table[i].step_valence)
            // Avalon❗✔️:             {
            let mut value = entry.from_valence;
            while value <= entry.to_valence {
                // Avalon❗✔️:                h=val-bond_electrons+charge;
                let hydrogens = value - bond_electrons + charge;
                // Avalon❗✔️:                if (h < 0) continue;
                if hydrogens < 0 {
                    value += entry.step_valence;
                    continue;
                }
                // Avalon❗✔️:                if (valence_table[i].lone_pairs > 0)
                // Avalon❗✔️:                   return (h);
                // Avalon❗✔️:                else
                // Avalon❗✔️:                   return (0);
                return if entry.lone_pairs > 0 { hydrogens } else { 0 };
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:          }
            // Avalon❗✔️:          else /* charge < 0 */
            // Avalon❗✔️:          {
        } else {
            // Avalon❗✔️:             for (val=valence_table[i].from_valence;
            // Avalon❗✔️:                  val<=valence_table[i].to_valence;
            // Avalon❗✔️:                  val+=valence_table[i].step_valence)
            // Avalon❗✔️:             {
            let mut value = entry.from_valence;
            while value <= entry.to_valence {
                // Avalon❗✔️:                h=val-bond_electrons-charge;
                let mut hydrogens = value - bond_electrons - charge;
                // Avalon❗✔️:                if (h < 0) continue;
                if hydrogens < 0 {
                    value += entry.step_valence;
                    continue;
                }
                // Avalon❗✔️:                if (valence_table[i].pair_deficit<h)
                // Avalon❗✔️:                       h=valence_table[i].pair_deficit;
                if entry.pair_deficit < hydrogens {
                    hydrogens = entry.pair_deficit;
                }
                // Avalon❗✔️:                /* hydrid ions are not stable if there are lone-pairs. */
                // Avalon❗✔️:                if (valence_table[i].lone_pairs > 0)
                // Avalon❗✔️:                   return (0);
                // Avalon❗✔️:                else
                // Avalon❗✔️:                   return (h);
                return if entry.lone_pairs > 0 { 0 } else { hydrogens };
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    return(0);
    // Avalon❗✔️: }
    0
}

fn compute_implicit_h(molecule: &MoleculeState, h_count: &mut [i32]) -> Result<(), FingerprintError> {
    // Avalon❗✔️: void
    // Avalon❗✔️: ComputeImplicitH(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                  int H_count[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Computes the implicit hydrogen counts for the atoms in *mp. The
    // Avalon❗✔️:  * result is returned in H_count[] starting at H_count[1] for atom
    // Avalon❗✔️:  * number 1. H_count[0] is not changed.
    // Avalon❗✔️:  *
    // Avalon❗✔️:  * Note: The elements of H_count are changed only if the are zero!!
    // Avalon❗✔️:  *       They need to be initialized by the caller!
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    if h_count.len() != molecule.atoms.len() + 1 {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon H-count array has the wrong one-based length".to_string(),
        });
    }
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    /* static */
    // Avalon❗✔️:    int *single_bond,      /* Counts of attached bonds of a type, */
    // Avalon❗✔️:        *aromatic_bond,    /* <array>[0] is unused */
    // Avalon❗✔️:        *double_bond,
    // Avalon❗✔️:        *triple_bond,
    // Avalon❗✔️:        *radical,	  /* total count of property */
    // Avalon❗✔️:        *charge;
    // Avalon❗✔️:    struct reaccs_atom_t *ap;
    // Avalon❗✔️:    struct reaccs_bond_t *bp;
    // Avalon❗✔️:
    let array_len = molecule.atoms.len() + 1;
    // Avalon❗✔️:    single_bond   = TypeAlloc(mp->n_atoms+1, int);
    let mut single_bond = vec![0_i32; array_len];
    // Avalon❗✔️:    aromatic_bond = TypeAlloc(mp->n_atoms+1, int);
    let mut aromatic_bond = vec![0_i32; array_len];
    // Avalon❗✔️:    double_bond   = TypeAlloc(mp->n_atoms+1, int);
    let mut double_bond = vec![0_i32; array_len];
    // Avalon❗✔️:    triple_bond   = TypeAlloc(mp->n_atoms+1, int);
    let mut triple_bond = vec![0_i32; array_len];
    // Avalon❗✔️:    radical       = TypeAlloc(mp->n_atoms+1, int);
    let mut radical = vec![0_i32; array_len];
    // Avalon❗✔️:    charge        = TypeAlloc(mp->n_atoms+1, int);
    let mut charge = vec![0_i32; array_len];
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:    {
    for (index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:       radical[i+1] = mp->atom_array[i].radical == DOUBLET;
        radical[index + 1] = i32::from(atom.radical == DOUBLET);
        // Avalon❗✔️:       charge[i+1]  = mp->atom_array[i].charge;
        charge[index + 1] = atom.charge;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<=mp->n_atoms; i++)
    // Avalon❗✔️:       single_bond[i] = aromatic_bond[i] = double_bond[i] = triple_bond[i] = 0;
    // Rust zero-initializes the four arrays at allocation.
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0,bp=mp->bond_array; i<mp->n_bonds; i++,bp++)
    // Avalon❗✔️:    {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        let at0 = source_atom_number(bond.atoms[0], molecule.atoms.len(), bond_index)?;
        let at1 = source_atom_number(bond.atoms[1], molecule.atoms.len(), bond_index)?;
        // Avalon❗✔️:       switch (bp->bond_type)
        // Avalon❗✔️:       {
        match bond.bond_type {
            // Avalon❗✔️:          case SINGLE: single_bond[bp->atoms[0]]++;
            // Avalon❗✔️:                       single_bond[bp->atoms[1]]++;
            // Avalon❗✔️:                       break;
            SINGLE => {
                single_bond[at0] += 1;
                single_bond[at1] += 1;
            }
            // Avalon❗✔️:          case DOUBLE: double_bond[bp->atoms[0]]++;
            // Avalon❗✔️:                       double_bond[bp->atoms[1]]++;
            // Avalon❗✔️:                       break;
            DOUBLE => {
                double_bond[at0] += 1;
                double_bond[at1] += 1;
            }
            // Avalon❗✔️:          case TRIPLE: triple_bond[bp->atoms[0]]++;
            // Avalon❗✔️:                       triple_bond[bp->atoms[1]]++;
            // Avalon❗✔️:                       break;
            TRIPLE => {
                triple_bond[at0] += 1;
                triple_bond[at1] += 1;
            }
            // Avalon❗✔️:          case AROMATIC: aromatic_bond[bp->atoms[0]]++;
            // Avalon❗✔️:                         aromatic_bond[bp->atoms[1]]++;
            // Avalon❗✔️:                         break;
            AROMATIC => {
                aromatic_bond[at0] += 1;
                aromatic_bond[at1] += 1;
            }
            // Avalon❗✔️:          default :
            // Avalon❗✔️:             single_bond[bp->atoms[0]]++;
            // Avalon❗✔️:             single_bond[bp->atoms[1]]++;
            // Avalon❗✔️:             break;
            _ => {
                single_bond[at0] += 1;
                single_bond[at1] += 1;
            } // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0,ap=mp->atom_array; i<mp->n_atoms; i++,ap++)
    // Avalon❗✔️:       if (H_count[i+1] == 0)
    // Avalon❗✔️:       {
    for (index, atom) in molecule.atoms.iter().enumerate() {
        if h_count[index + 1] == 0 {
            // Avalon❗✔️: 	 H_count[i+1] = ImplicitHydrogens(ap->atom_symbol,
            // Avalon❗✔️: 					  single_bond[i+1],
            // Avalon❗✔️: 					  aromatic_bond[i+1],
            // Avalon❗✔️: 					  double_bond[i+1],
            // Avalon❗✔️: 					  triple_bond[i+1],
            // Avalon❗✔️: 					  radical[i+1],
            // Avalon❗✔️: 					  charge[i+1]);
            h_count[index + 1] = implicit_hydrogens(
                &atom.atom_symbol,
                single_bond[index + 1],
                aromatic_bond[index + 1],
                double_bond[index + 1],
                triple_bond[index + 1],
                radical[index + 1],
                charge[index + 1],
            );
            // Avalon❗✔️: 	 if (H_count[i+1] < 0) H_count[i+1] = 0;
            if h_count[index + 1] < 0 {
                h_count[index + 1] = 0;
            }
            // Avalon❗✔️:       }
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    MyFree((char *)single_bond); MyFree((char *)aromatic_bond);
    // Avalon❗✔️:    MyFree((char *)double_bond); MyFree((char *)triple_bond);
    // Avalon❗✔️:    MyFree((char *)radical); MyFree((char *)charge);
    // Rust vectors release the same six working arrays on scope exit.
    // Avalon❗✔️: }
    Ok(())
}

fn source_atom_number(source_index: i32, atom_count: usize, bond_index: usize) -> Result<usize, FingerprintError> {
    usize::try_from(source_index)
        .ok()
        .filter(|&value| value > 0 && value <= atom_count)
        .ok_or_else(|| FingerprintError::AvalonConversion {
            reason: format!(
                "Avalon bond {} references invalid atom {}",
                bond_index + 1,
                source_index
            ),
        })
}

fn guess_h_counts_from_substitution(
    molecule: &mut MoleculeState,
    neighbours: &[Neighbourhood],
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: void GuessHCountsFromSubstitution(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                                   neighbourhood_t nbp[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Uses substitution count query options to guess the number of
    // Avalon❗✔️:  * required hydrogens on oxigen and nitrogen atoms. This can speed
    // Avalon❗✔️:  * up screening for substructure matches.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i, j, nsingle, ndouble, ntriple, naromatic, nother, nexplicit;
    // Avalon❗✔️:    struct reaccs_atom_t *ap;
    // Avalon❗✔️:
    // Avalon❗✔️:    if (!mp) return;
    if neighbours.len() != molecule.atoms.len() {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon neighbourhood array has the wrong length".to_string(),
        });
    }
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for (index, neighbourhood) in neighbours.iter().enumerate() {
        let atom = &molecule.atoms[index];
        // Avalon❗✔️:       /* ignore if already defined */
        // Avalon❗✔️:       if (ap->query_H_count != NONE) continue;
        if atom.query_h_count != NONE {
            continue;
        }
        // Avalon❗✔️:       /* Only easy cases */
        // Avalon❗✔️:       if (ap->charge != 0) continue;
        if atom.charge != 0 {
            continue;
        }
        // Avalon❗✔️:       if (ap->radical != 0) continue;
        if atom.radical != 0 {
            continue;
        }
        // Avalon❗✔️:       nsingle = ndouble = naromatic= ntriple = nother = 0;
        let (mut nsingle, mut ndouble, mut naromatic, mut ntriple, mut nother) = (0_i32, 0_i32, 0_i32, 0_i32, 0_i32);
        // Avalon❗✔️:       nexplicit = 0;
        let mut nexplicit = 0_i32;
        // Avalon❗✔️:       for (j=0; j<nbp[i].n_ligands; j++)
        // Avalon❗✔️:          if (AtomSymbolMatch(mp->atom_array[nbp[i].atoms[j]].atom_symbol, "H,D,T"))
        // Avalon❗✔️:             nexplicit++;
        for &neighbour_atom in neighbourhood.atoms() {
            if matches!(molecule.atoms[neighbour_atom].atom_symbol.as_str(), "H" | "D" | "T") {
                nexplicit += 1;
            }
        }
        // Avalon❗✔️:       /* Collect bond orders of ligands */
        // Avalon❗✔️:       for (j=0; j<nbp[i].n_ligands; j++)
        // Avalon❗✔️:          switch (mp->bond_array[nbp[i].bonds[j]].bond_type)
        // Avalon❗✔️:          {
        for &bond_index in neighbourhood.bonds() {
            match molecule.bonds[bond_index].bond_type {
                // Avalon❗✔️:             case SINGLE: nsingle++; break;
                SINGLE => nsingle += 1,
                // Avalon❗✔️:             case DOUBLE: ndouble++; break;
                DOUBLE => ndouble += 1,
                // Avalon❗✔️:             case AROMATIC: naromatic++; break;
                AROMATIC => naromatic += 1,
                // Avalon❗✔️:             case TRIPLE: ntriple++; break;
                TRIPLE => ntriple += 1,
                // Avalon❗✔️:             default: nother++;
                _ => nother += 1,
                // Avalon❗✔️:          }
            }
        }
        // Avalon❗✔️:       /* Ignore if there are query bonds */
        // Avalon❗✔️:       if (nother > 0) continue;
        if nother > 0 {
            continue;
        }
        let atom = &mut molecule.atoms[index];
        // Avalon❗✔️:       /* for now, only handle AS_IS case */
        // Avalon❗✔️:       if (0 == strcmp("C", ap->atom_symbol))
        // Avalon❗✔️:       {
        if atom.atom_symbol == "C" {
            // Avalon❗✔️:          if (naromatic == 2)
            // Avalon❗✔️:          {
            if naromatic == 2 {
                // Avalon❗✔️:             nsingle++; ndouble++; naromatic=0;
                nsingle += 1;
                ndouble += 1;
                naromatic = 0;
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          if (naromatic > 0) continue;
            if naromatic > 0 {
                continue;
            }
            // Avalon❗✔️:          if (ap->sub_desc == NONE)      /* degree from explicit hydrogens */
            // Avalon❗✔️:          {
            if atom.sub_desc == NONE {
                // Avalon❗✔️:             if (nsingle+2*ndouble+3*ntriple == 4)
                // Avalon❗✔️:             {
                if nsingle + 2 * ndouble + 3 * ntriple == 4 {
                    // Avalon❗✔️: // fprintf(stderr, "%d: s=%d, d=%d, t=%d, sum=%d\n",
                    // Avalon❗✔️: // i+1, nsingle, ndouble, ntriple, nsingle+2*ndouble+3*ntriple);
                    // Avalon❗✔️:                ap->sub_desc = SUB_AS_IS;
                    atom.sub_desc = SUB_AS_IS;
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
                // Avalon❗✔️:          else if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
                // Avalon❗✔️:          {
            } else if atom.sub_desc == SUB_AS_IS {
                // Avalon❗✔️:             ap->query_H_count =
                // Avalon❗✔️:                nexplicit+ZERO_COUNT+4-nsingle-2*ndouble-3*ntriple;
                atom.query_h_count = nexplicit + ZERO_COUNT + 4 - nsingle - 2 * ndouble - 3 * ntriple;
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (0 == strcmp("O", ap->atom_symbol))
            // Avalon❗✔️:       {
        } else if atom.atom_symbol == "O" {
            // Avalon❗✔️:          if (ntriple+naromatic != 0) continue;
            if ntriple + naromatic != 0 {
                continue;
            }
            // Avalon❗✔️:          if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
            // Avalon❗✔️:             ap->query_H_count = nexplicit+ZERO_COUNT+2-nsingle-2*ndouble;
            if atom.sub_desc == SUB_AS_IS {
                atom.query_h_count = nexplicit + ZERO_COUNT + 2 - nsingle - 2 * ndouble;
            }
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (0 == strcmp("N", ap->atom_symbol))
            // Avalon❗✔️:       {
        } else if atom.atom_symbol == "N" {
            // Avalon❗✔️:          if (naromatic != 0) continue;
            if naromatic != 0 {
                continue;
            }
            // Avalon❗✔️:          if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
            // Avalon❗✔️:             ap->query_H_count = nexplicit+ZERO_COUNT+3-nsingle-2*ndouble;
            if atom.sub_desc == SUB_AS_IS {
                atom.query_h_count = nexplicit + ZERO_COUNT + 3 - nsingle - 2 * ndouble;
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️: // fprintf(stderr,"UTILITIES:\t'%s' atom %d has h_count=%d\n", ap->atom_symbol, i+1, ap->query_H_count);
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️: }
    Ok(())
}

pub(super) fn collect_hydrogen_counts(
    molecule: &mut MoleculeState,
    neighbours: &[Neighbourhood],
    as_query: bool,
) -> Result<Vec<i32>, FingerprintError> {
    // Avalon❗✔️:    /* collect hydrogen counts */
    // Avalon❗✔️:    H_count = TypeAlloc(mp->n_atoms+1, int);
    let mut h_count = vec![0_i32; molecule.atoms.len() + 1];
    // Avalon❗✔️:    if (as_query)
    // Avalon❗✔️:    {
    if as_query {
        // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        // Avalon❗✔️:       {
        for (bond_index, bond) in molecule.bonds.iter().enumerate() {
            let at0 = source_atom_number(bond.atoms[0], molecule.atoms.len(), bond_index)?;
            let at1 = source_atom_number(bond.atoms[1], molecule.atoms.len(), bond_index)?;
            // Avalon❗✔️:          if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            // Avalon❗✔️:             H_count[bp->atoms[1]]++;
            if molecule.atoms[at0 - 1].atom_symbol == "H" {
                h_count[at1] += 1;
                // Avalon❗✔️:          else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
                // Avalon❗✔️:             H_count[bp->atoms[0]]++;
            } else if molecule.atoms[at1 - 1].atom_symbol == "H" {
                h_count[at0] += 1;
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       GuessHCountsFromSubstitution(mp, nbp);
        guess_h_counts_from_substitution(molecule, neighbours)?;
        // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
        // Avalon❗✔️:          if (ap->query_H_count != NONE)
        // Avalon❗✔️:          {
        for (index, atom) in molecule.atoms.iter().enumerate() {
            if atom.query_h_count != NONE {
                // Avalon❗✔️:             H_count[i+1] = ap->query_H_count-ZERO_COUNT;
                h_count[index + 1] = atom.query_h_count - ZERO_COUNT;
                // Avalon❗✔️: // fprintf(stderr,"SSMATCH:\t'%s' atom %d has h_count=%d\n", ap->atom_symbol, i+1, ap->query_H_count);
                // Avalon❗✔️:          }
            }
        }
        // Avalon❗✔️:    }
        // Avalon❗✔️:    else
        // Avalon❗✔️:    {
    } else {
        // Avalon❗✔️:       ComputeImplicitH(mp, H_count);
        compute_implicit_h(molecule, &mut h_count)?;
        // Avalon❗✔️:       /* Add the explicit hydrogens to the implicit counts */
        // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        // Avalon❗✔️:       {
        for (bond_index, bond) in molecule.bonds.iter().enumerate() {
            let at0 = source_atom_number(bond.atoms[0], molecule.atoms.len(), bond_index)?;
            let at1 = source_atom_number(bond.atoms[1], molecule.atoms.len(), bond_index)?;
            // Avalon❗✔️:          if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            // Avalon❗✔️:             H_count[bp->atoms[1]]++;
            if molecule.atoms[at0 - 1].atom_symbol == "H" {
                h_count[at1] += 1;
                // Avalon❗✔️:          else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
                // Avalon❗✔️:             H_count[bp->atoms[0]]++;
            } else if molecule.atoms[at1 - 1].atom_symbol == "H" {
                h_count[at0] += 1;
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    Ok(h_count)
}

pub(super) fn ring_state(molecule: &MoleculeState) -> Result<(Vec<i32>, Vec<i32>), FingerprintError> {
    // Avalon❗✔️: typedef unsigned atom_pair[2];
    // Avalon❗✔️: void RingState(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                int atom_status[],
    // Avalon❗✔️:                int bond_status[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Computes how many basis rings each bond shares and how many
    // Avalon❗✔️:  * ring bonds are attached to an atom. The results are stored in
    // Avalon❗✔️:  * atom_status[] and bond_status[] respectively.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    bond_set_node *rph, *ring_list;
    // Avalon❗✔️:    atom_pair *bonds;
    // Avalon❗✔️:    struct reaccs_bond_t *bp;
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:       atom_status[i] = 0;
    let mut atom_status = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:       bond_status[i] = 0;
    let mut bond_status = vec![0_i32; molecule.bonds.len()];
    // Avalon❗✔️:
    // Avalon❗✔️:    if (mp->n_bonds == 0) return;
    if molecule.bonds.is_empty() {
        return Ok((atom_status, bond_status));
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    bonds = TypeAlloc(mp->n_bonds, atom_pair); /* get basis rings */
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:    {
    let mut bonds = Vec::with_capacity(molecule.bonds.len());
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Validate the same one-based source indices before entering the
        // source graph routine, whose out-of-range behavior is undefined.
        source_atom_number(bond.atoms[0], molecule.atoms.len(), bond_index)?;
        source_atom_number(bond.atoms[1], molecule.atoms.len(), bond_index)?;
        // Avalon❗✔️:       bonds[i][0] = mp->bond_array[i].atoms[0];
        // Avalon❗✔️:       bonds[i][1] = mp->bond_array[i].atoms[1];
        bonds.push([bond.atoms[0] as usize, bond.atoms[1] as usize]);
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    ring_list = RingList(bonds,mp->n_bonds);
    let mut rings = ring_list(&bonds);
    // Avalon❗✔️:    ring_list = CombineRings(ring_list);
    combine_rings(&mut rings);
    // Avalon❗✔️:    MyFree((char *)bonds);
    // Rust releases the temporary bond pairs on scope exit.
    // Avalon❗✔️:
    // Avalon❗✔️:    for (rph=ring_list; rph; rph=rph->next)
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:          if (IsMember(rph->bond_set,i))
    // Avalon❗✔️:           bond_status[i]++;
    for ring in &rings {
        for (bond_index, status) in bond_status.iter_mut().enumerate() {
            if ring.bond_set.contains(bond_index) {
                *status += 1;
            }
        }
    }
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       if (bond_status[i] > 0)
    // Avalon❗✔️:       {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        if bond_status[bond_index] > 0 {
            // Avalon❗✔️:          atom_status[bp->atoms[0]-1]++;
            atom_status[bond.atoms[0] as usize - 1] += 1;
            // Avalon❗✔️:          atom_status[bp->atoms[1]-1]++;
            atom_status[bond.atoms[1] as usize - 1] += 1;
            // Avalon❗✔️:       }
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    DisposeBondSetList(ring_list);
    // Rust releases the basis-ring vector and its sets on scope exit.
    // Avalon❗✔️: }
    Ok((atom_status, bond_status))
}

#[allow(clippy::too_many_arguments)]
fn mark_recursive(
    molecule: &mut MoleculeState,
    touched_atoms: &mut [bool],
    touched_bonds: &mut [bool],
    start_index: usize,
    path_length: usize,
    current_index: usize,
    max_size: usize,
    neighbours: &[Neighbourhood],
) {
    // Avalon❗✔️: static
    // Avalon❗✔️: void MarkRecursive(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                    int touched_atoms[],
    // Avalon❗✔️:                    int touched_bonds[],
    // Avalon❗✔️:                    int start_index,
    // Avalon❗✔️:                    int path_length,
    // Avalon❗✔️:                    int current_index,
    // Avalon❗✔️:                    int max_size,
    // Avalon❗✔️:                    neighbourhood_t nbp[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Recursively enumerates the rings in *mp and sets the rsize_flags bits
    // Avalon❗✔️:  * when a ring is found.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i, j, ai, bi;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<nbp[current_index].n_ligands; i++)
    // Avalon❗✔️:    {
    for ligand in 0..neighbours[current_index].n_ligands {
        // Avalon❗✔️:       ai = nbp[current_index].atoms[i];
        let atom_index = neighbours[current_index].atoms[ligand];
        // Avalon❗✔️:       /* rings only count if started at least atom index */
        // Avalon❗✔️:       if (ai < start_index) continue;
        if atom_index < start_index {
            continue;
        }
        // Avalon❗✔️:       if (ai == start_index)                    /* ring found */
        // Avalon❗✔️:       {
        if atom_index == start_index {
            // Avalon❗✔️:          if (path_length < 3) continue;
            if path_length < 3 {
                continue;
            }
            // Avalon❗✔️:          for (j=0; j<mp->n_atoms; j++)
            // Avalon❗✔️:             if (touched_atoms[j])
            // Avalon❗✔️:                mp->atom_array[j].rsize_flags |= (1<<path_length);
            for (index, &touched) in touched_atoms.iter().enumerate() {
                if touched {
                    molecule.atoms[index].rsize_flags |= 1_u32 << path_length;
                }
            }
            // Avalon❗✔️:          for (j=0; j<mp->n_bonds; j++)
            // Avalon❗✔️:             if (touched_bonds[j])
            // Avalon❗✔️:                mp->bond_array[j].rsize_flags |= (1<<path_length);
            for (index, &touched) in touched_bonds.iter().enumerate() {
                if touched {
                    molecule.bonds[index].rsize_flags |= 1_u32 << path_length;
                }
            }
            // Avalon❗✔️:          continue;
            continue;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (touched_atoms[ai]) continue;   /* don't walk backwards */
        if touched_atoms[atom_index] {
            continue;
        }
        // Avalon❗✔️:       /* continue recursion */
        // Avalon❗✔️:       if (path_length+1 > max_size)             /* don't go too far */
        // Avalon❗✔️:          continue;
        if path_length + 1 > max_size {
            continue;
        }
        // Avalon❗✔️:       if (mp->atom_array[ai].rsize_flags == 0)  /* only ring atoms */
        // Avalon❗✔️:          continue;
        if molecule.atoms[atom_index].rsize_flags == 0 {
            continue;
        }
        // Avalon❗✔️:       bi = nbp[current_index].bonds[i];
        let bond_index = neighbours[current_index].bonds[ligand];
        // Avalon❗✔️:       if (mp->bond_array[bi].rsize_flags == 0)  /* only ring bonds */
        // Avalon❗✔️:          continue;
        if molecule.bonds[bond_index].rsize_flags == 0 {
            continue;
        }
        // Avalon❗✔️:       touched_atoms[ai] = 1;  /* updating */
        touched_atoms[atom_index] = true;
        // Avalon❗✔️:       touched_bonds[bi] = 1;
        touched_bonds[bond_index] = true;
        // Avalon❗✔️:       MarkRecursive(mp, touched_atoms, touched_bonds,
        // Avalon❗✔️:                     start_index, path_length+1, ai, max_size, nbp);
        mark_recursive(
            molecule,
            touched_atoms,
            touched_bonds,
            start_index,
            path_length + 1,
            atom_index,
            max_size,
            neighbours,
        );
        // Avalon❗✔️:       // touched_bonds[nbp[current_index].bonds[i]] = 0;
        // Avalon❗✔️:       touched_atoms[ai] = 0;  /* down-dating */
        touched_atoms[atom_index] = false;
        // Avalon❗✔️:       touched_bonds[bi] = 0;
        touched_bonds[bond_index] = false;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️: }
}

pub(super) fn set_ring_size_flags(
    molecule: &mut MoleculeState,
    max_size: usize,
    neighbours: &[Neighbourhood],
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: void SetRingSizeFlags(struct reaccs_molecule_t *mp, int max_size,
    // Avalon❗✔️:                       neighbourhood_t nbp[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Sets the perceived ring-sizes as flag bits in the bond rsize_flags field.
    // Avalon❗✔️:  * It used a recursive enumeration algorithm pruned to only ring bonds/atoms.
    // Avalon❗✔️:  * nbp[] is the neighbourhood array to speed up processing.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int *atom_status, *bond_status;
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    int *touched_atoms; /* this array is up- and down-dated during recursion */
    // Avalon❗✔️:    int *touched_bonds; /* this array is up- and down-dated during recursion */
    if neighbours.len() != molecule.atoms.len() {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon neighbourhood array has the wrong length".to_string(),
        });
    }
    if max_size >= u32::BITS as usize {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon ring-size bit limit exceeds the source field width".to_string(),
        });
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    /* label cyclic parts with unused rsize_flags bit */
    // Avalon❗✔️:    atom_status = TypeAlloc(mp->n_atoms, int);
    // Avalon❗✔️:    bond_status = TypeAlloc(mp->n_bonds, int);
    // Avalon❗✔️:    RingState(mp, atom_status, bond_status);
    let (atom_status, bond_status) = ring_state(molecule)?;
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:       if (atom_status[i] > 0) mp->atom_array[i].rsize_flags = 1;
    // Avalon❗✔️:       else                    mp->atom_array[i].rsize_flags = 0;
    for (atom, status) in molecule.atoms.iter_mut().zip(atom_status) {
        atom.rsize_flags = u32::from(status > 0);
    }
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:       if (bond_status[i] > 0) mp->bond_array[i].rsize_flags = 1;
    // Avalon❗✔️:       else                    mp->bond_array[i].rsize_flags = 0;
    for (bond, status) in molecule.bonds.iter_mut().zip(bond_status) {
        bond.rsize_flags = u32::from(status > 0);
    }
    // Avalon❗✔️:    if (atom_status) MyFree((char *)atom_status);
    // Avalon❗✔️:    if (bond_status) MyFree((char *)bond_status);
    // Rust releases both status vectors after their zip loops.
    // Avalon❗✔️:
    // Avalon❗✔️:    touched_atoms = TypeAlloc(mp->n_atoms, int);
    let mut touched_atoms = vec![false; molecule.atoms.len()];
    // Avalon❗✔️:    touched_bonds = TypeAlloc(mp->n_bonds, int);
    let mut touched_bonds = vec![false; molecule.bonds.len()];
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:    {
    for index in 0..molecule.atoms.len() {
        // Avalon❗✔️:       if (mp->atom_array[i].rsize_flags == 0) continue;
        if molecule.atoms[index].rsize_flags == 0 {
            continue;
        }
        // Avalon❗✔️:       touched_atoms[i]=1;     /* updating */
        touched_atoms[index] = true;
        // Avalon❗✔️:       MarkRecursive(mp, touched_atoms, touched_bonds, i, 1, i, max_size, nbp);
        mark_recursive(
            molecule,
            &mut touched_atoms,
            &mut touched_bonds,
            index,
            1,
            index,
            max_size,
            neighbours,
        );
        // Avalon❗✔️:       touched_atoms[i]=0;     /* down-dating */
        touched_atoms[index] = false;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    MyFree((char *)touched_atoms);
    // Avalon❗✔️:    MyFree((char *)touched_bonds);
    // Rust releases both traversal vectors on scope exit.
    // Avalon❗✔️: }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    fn state_with_bonds(atom_count: usize, endpoints: &[[i32; 2]]) -> MoleculeState {
        MoleculeState {
            atoms: vec![Atom::default(); atom_count],
            bonds: endpoints
                .iter()
                .map(|&atoms| Bond {
                    atoms,
                    ..Bond::default()
                })
                .collect(),
            ..MoleculeState::default()
        }
    }

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

    #[test]
    fn setup_neighbourhood_preserves_source_bond_table_order() {
        let molecule = state_with_bonds(4, &[[1, 3], [2, 1], [1, 4], [3, 2]]);
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        assert_eq!(neighbours[0].atoms(), &[2, 1, 3]);
        assert_eq!(neighbours[0].bonds(), &[0, 1, 2]);
        assert_eq!(neighbours[1].atoms(), &[0, 2]);
        assert_eq!(neighbours[1].bonds(), &[1, 3]);
        assert_eq!(neighbours[2].atoms(), &[0, 1]);
        assert_eq!(neighbours[2].bonds(), &[0, 3]);
        assert_eq!(neighbours[3].atoms(), &[0]);
        assert_eq!(neighbours[3].bonds(), &[2]);
    }

    #[test]
    fn setup_neighbourhood_applies_source_nlimit_to_both_endpoints() {
        let molecule = state_with_bonds(4, &[[1, 2], [2, 3], [1, 4]]);
        let neighbours = setup_neighbourhood(&molecule, 2).unwrap();

        assert_eq!(neighbours[0].atoms(), &[1]);
        assert_eq!(neighbours[1].atoms(), &[0]);
        assert!(neighbours[2].atoms().is_empty());
        assert!(neighbours[3].atoms().is_empty());
    }

    #[test]
    fn setup_neighbourhood_accepts_twenty_and_rejects_twenty_one_ligands() {
        let mut endpoints = Vec::new();
        for atom in 2..=21 {
            endpoints.push([1, atom]);
        }
        let molecule = state_with_bonds(21, &endpoints);
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();
        assert_eq!(neighbours[0].atoms().len(), MAX_NEIGHBOURS);

        endpoints.push([1, 22]);
        let molecule = state_with_bonds(22, &endpoints);
        assert!(matches!(
            setup_neighbourhood(&molecule, molecule.atoms.len()),
            Err(FingerprintError::AvalonConversion { reason })
                if reason == "Avalon neighbourhood capacity exceeded at atom 1"
        ));
    }

    #[test]
    fn setup_neighbourhood_rejects_invalid_one_based_atom_indices() {
        for invalid in [0, -1, 3] {
            let molecule = state_with_bonds(2, &[[1, invalid]]);
            assert!(matches!(
                setup_neighbourhood(&molecule, molecule.atoms.len()),
                Err(FingerprintError::AvalonConversion { .. })
            ));
        }
    }

    #[test]
    fn implicit_hydrogens_uses_source_valence_and_aromatic_electron_rules() {
        assert_eq!(implicit_hydrogens("C", 0, 0, 0, 0, 0, 0), 4);
        assert_eq!(implicit_hydrogens("C", 1, 0, 0, 0, 0, 0), 3);
        assert_eq!(implicit_hydrogens("N", 0, 0, 0, 0, 0, 0), 3);
        assert_eq!(implicit_hydrogens("N", 0, 2, 0, 0, 0, 0), 0);
        assert_eq!(implicit_hydrogens("C", 0, 4, 0, 0, 0, 0), 0);
        assert_eq!(implicit_hydrogens("Xe", 0, 0, 0, 0, 0, 0), 0);
    }

    #[test]
    fn compute_implicit_h_preserves_nonzero_one_based_slots() {
        let molecule = MoleculeState {
            atoms: vec![atom("C"), atom("C")],
            bonds: vec![bond([1, 2], SINGLE)],
            ..MoleculeState::default()
        };
        let mut h_count = vec![77, 2, 0];

        compute_implicit_h(&molecule, &mut h_count).unwrap();

        assert_eq!(h_count, vec![77, 2, 3]);
    }

    #[test]
    fn ordinary_counts_add_explicit_h_after_implicit_h() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"), atom("H")],
            bonds: vec![bond([1, 2], SINGLE)],
            ..MoleculeState::default()
        };
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        let h_count = collect_hydrogen_counts(&mut molecule, &neighbours, false).unwrap();

        assert_eq!(h_count, vec![0, 4, 0]);
    }

    #[test]
    fn query_h_count_overrides_precounted_explicit_h() {
        let mut query_atom = atom("N");
        query_atom.query_h_count = ZERO_COUNT + 3;
        let mut molecule = MoleculeState {
            atoms: vec![query_atom, atom("H")],
            bonds: vec![bond([1, 2], SINGLE)],
            ..MoleculeState::default()
        };
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        let h_count = collect_hydrogen_counts(&mut molecule, &neighbours, true).unwrap();

        assert_eq!(h_count, vec![0, 3, 0]);
    }

    #[test]
    fn substitution_guess_counts_hydrogen_isotopes_but_rejects_query_bonds() {
        let mut oxygen = atom("O");
        oxygen.sub_desc = SUB_AS_IS;
        let mut molecule = MoleculeState {
            atoms: vec![oxygen, atom("C"), atom("D")],
            bonds: vec![bond([1, 2], SINGLE), bond([1, 3], SINGLE)],
            ..MoleculeState::default()
        };
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();
        guess_h_counts_from_substitution(&mut molecule, &neighbours).unwrap();
        assert_eq!(molecule.atoms[0].query_h_count, ZERO_COUNT + 1);

        molecule.atoms[0].query_h_count = NONE;
        molecule.bonds[1].bond_type = 5;
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();
        guess_h_counts_from_substitution(&mut molecule, &neighbours).unwrap();
        assert_eq!(molecule.atoms[0].query_h_count, NONE);
    }

    #[test]
    fn carbon_none_substitution_only_transitions_to_as_is_on_first_pass() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"), atom("C"), atom("C"), atom("C"), atom("C")],
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([1, 3], SINGLE),
                bond([1, 4], SINGLE),
                bond([1, 5], SINGLE),
            ],
            ..MoleculeState::default()
        };
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        guess_h_counts_from_substitution(&mut molecule, &neighbours).unwrap();

        assert_eq!(molecule.atoms[0].sub_desc, SUB_AS_IS);
        assert_eq!(molecule.atoms[0].query_h_count, NONE);
    }

    #[test]
    fn ring_state_counts_basis_membership_and_attached_ring_bonds() {
        let molecule = MoleculeState {
            atoms: vec![atom("C"); 4],
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([2, 3], SINGLE),
                bond([3, 4], SINGLE),
                bond([4, 1], SINGLE),
                bond([1, 3], SINGLE),
            ],
            ..MoleculeState::default()
        };

        let (atom_status, bond_status) = ring_state(&molecule).unwrap();

        assert_eq!(atom_status, vec![3, 2, 3, 2]);
        assert_eq!(bond_status, vec![1, 1, 1, 1, 2]);
    }

    #[test]
    fn ring_size_flags_enumerate_triangles_and_outer_square() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 4],
            bonds: vec![
                bond([1, 2], SINGLE),
                bond([2, 3], SINGLE),
                bond([3, 4], SINGLE),
                bond([4, 1], SINGLE),
                bond([1, 3], SINGLE),
            ],
            ..MoleculeState::default()
        };
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        set_ring_size_flags(&mut molecule, 14, &neighbours).unwrap();

        let triangle_and_square = 1 | (1 << 3) | (1 << 4);
        assert!(
            molecule
                .atoms
                .iter()
                .all(|atom| atom.rsize_flags == triangle_and_square)
        );
        assert_eq!(
            molecule.bonds.iter().map(|bond| bond.rsize_flags).collect::<Vec<_>>(),
            vec![
                triangle_and_square,
                triangle_and_square,
                triangle_and_square,
                triangle_and_square,
                1 | (1 << 3),
            ]
        );
    }

    #[test]
    fn ring_size_flags_clear_stale_flags_on_acyclic_graph() {
        let mut molecule = MoleculeState {
            atoms: vec![atom("C"); 3],
            bonds: vec![bond([1, 2], SINGLE), bond([2, 3], SINGLE)],
            ..MoleculeState::default()
        };
        for atom in &mut molecule.atoms {
            atom.rsize_flags = u32::MAX;
        }
        for bond in &mut molecule.bonds {
            bond.rsize_flags = u32::MAX;
        }
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();

        set_ring_size_flags(&mut molecule, 14, &neighbours).unwrap();

        assert!(molecule.atoms.iter().all(|atom| atom.rsize_flags == 0));
        assert!(molecule.bonds.iter().all(|bond| bond.rsize_flags == 0));
    }
}
