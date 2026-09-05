//! Avalon atom-symbol predicates and source-order symbol-list lookup.

use super::reaccs::{MoleculeState, SymbolList};

// Avalon❗✔️: char *HC_table[] =                      /* pseudosymbol "G" */
// Avalon❗✔️:    {"H", "C", (char *)NULL};
const HC_TABLE: &[&str] = &["H", "C"];

// Avalon❗✔️: char *non_metal_hetero_elements[] =     /* pseudosymbol "Q" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "He",
// Avalon❗✔️:       "B", "N", "O", "F", "Ne",
// Avalon❗✔️:       "Si", "P", "S", "Cl", "Ar",
// Avalon❗✔️:       "As", "Se", "Br", "Kr",
// Avalon❗✔️:       "Sb", "Te", "I", "Xe",
// Avalon❗✔️:       "At",                     /* "Rn", This element must be removed */
// Avalon❗✔️:       (char *)NULL,             /* because of a trick in utils.c */
// Avalon❗✔️:    };
const NON_METAL_HETERO_ELEMENTS: &[&str] = &[
    "He", "B", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Sb", "Te", "I", "Xe", "At",
];

// Avalon❗✔️: char *metals[] =                /* pseudosymbol "M" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Li", "Be",
// Avalon❗✔️:       "Na", "Mg", "Al",
// Avalon❗✔️:       "K", "Ca", "Sc",
// Avalon❗✔️:         "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
// Avalon❗✔️:            "Ga",
// Avalon❗✔️:       "Rb", "Sr", "Y",
// Avalon❗✔️:         "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
// Avalon❗✔️:            "In", "Sn",
// Avalon❗✔️:       "Cs", "Ba", "La",
// Avalon❗✔️:        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
// Avalon❗✔️:        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
// Avalon❗✔️:         "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
// Avalon❗✔️:            "Tl", "Pb", "Bi", "Po",
// Avalon❗✔️:       "Fr", "Ra", "Ac",
// Avalon❗✔️:        "Th", "Pa", "U", "Np", "Pu",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const METALS: &[&str] = &[
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
];

// Avalon❗✔️: char *non_metal_small_solution[] =      /* pseudosymbol "Qs" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "H",
// Avalon❗✔️:       "B",  "C", "N", "O", "F",
// Avalon❗✔️:            "Si", "P", "S", "Cl",
// Avalon❗✔️:                      "Se", "Br",
// Avalon❗✔️:                             "I",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const NON_METAL_SMALL_SOLUTION: &[&str] = &["H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"];

// Avalon❗✔️: char *alkali_metals[] =                /* pseudosymbol "alk" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Li", "Na", "K", "Rb", "Cs", "Fr",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const ALKALI_METALS: &[&str] = &["Li", "Na", "K", "Rb", "Cs", "Fr"];

// Avalon❗✔️: char *gr2[] =                              /* pseudosymbol "gr2" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const GR2: &[&str] = &["Be", "Mg", "Ca", "Sr", "Ba", "Ra"];

// Avalon❗✔️: char *gr3[] =                               /* pseudosymbol "gr3" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "B", "Al", "Ga", "In", "Tl",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const GR3: &[&str] = &["B", "Al", "Ga", "In", "Tl"];

// Avalon❗✔️: char *gr4[] =                              /* pseudosymbol "gr4" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "C", "Si", "Ge", "Sn", "Pb",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const GR4: &[&str] = &["C", "Si", "Ge", "Sn", "Pb"];

// Avalon❗✔️: char *ONS_table[] =                     /* pseudosymbol "ONS" or "ons" */
// Avalon❗✔️:    {"O", "N", "S", (char *)NULL};
const ONS_TABLE: &[&str] = &["O", "N", "S"];

// Avalon❗✔️: char *on2[] =                            /* pseudosymbol "on2" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "O", "N", "S", "P", "Se", "Te", "Po",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:   };
const ON2: &[&str] = &["O", "N", "S", "P", "Se", "Te", "Po"];

// Avalon❗✔️: char *halogenes[] =                     /* pseudosymbol "X" or "hal" */
// Avalon❗✔️:    {"F", "Cl", "Br", "I", "At", (char *)NULL};
const HALOGENES: &[&str] = &["F", "Cl", "Br", "I", "At"];

// Avalon❗✔️: char *ha2[] =                     /* pseudosymbol "ha2" */
// Avalon❗✔️:    {"Cl", "Br", "I", "At", (char *)NULL};
const HA2: &[&str] = &["Cl", "Br", "I", "At"];

// Avalon❗✔️: char *transition_metals[] =                /* pseudosymbol "trn" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
// Avalon❗✔️:       "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
// Avalon❗✔️:       "La", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TRANSITION_METALS: &[&str] = &[
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
];

// Avalon❗✔️: char *tra[] =                /* pseudosymbol "tra" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
// Avalon❗✔️:       "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
// Avalon❗✔️:       "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TRA: &[&str] = &[
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt",
];

// Avalon❗✔️: char *trb[] =                /* pseudosymbol "trb" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Cu", "Zn", "Ag", "Cd", "Au", "Hg",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TRB: &[&str] = &["Cu", "Zn", "Ag", "Cd", "Au", "Hg"];

// Avalon❗✔️: char *tm1[] =                /* pseudosymbol "tm1" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Cu", "Ag", "Au",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM1: &[&str] = &["Cu", "Ag", "Au"];

// Avalon❗✔️: char *tm2[] =                /* pseudosymbol "tm2" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Zn", "Cd", "Hg",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM2: &[&str] = &["Zn", "Cd", "Hg"];

// Avalon❗✔️: char *tm3[] =                /* pseudosymbol "tm3" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Sc", "Y", "La",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM3: &[&str] = &["Sc", "Y", "La"];

// Avalon❗✔️: char *tm4[] =                /* pseudosymbol "tm4" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Ti", "Zr", "Hf",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM4: &[&str] = &["Ti", "Zr", "Hf"];

// Avalon❗✔️: char *tm5[] =                /* pseudosymbol "tm5" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "V", "Nb", "Ta",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM5: &[&str] = &["V", "Nb", "Ta"];

// Avalon❗✔️: char *tm6[] =                /* pseudosymbol "tm6" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Cr", "Mo", "W",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM6: &[&str] = &["Cr", "Mo", "W"];

// Avalon❗✔️: char *tm7[] =                /* pseudosymbol "tm7" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Mn", "Tc", "Re",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM7: &[&str] = &["Mn", "Tc", "Re"];

// Avalon❗✔️: char *tm8[] =                /* pseudosymbol "tm8" */
// Avalon❗✔️:    {
// Avalon❗✔️:       "Fe", "Co", "Ni",
// Avalon❗✔️:       "Ru", "Rh", "Pd",
// Avalon❗✔️:       "Os", "Ir", "Pt",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const TM8: &[&str] = &["Fe", "Co", "Ni", "Ru", "Rh", "Pd", "Os", "Ir", "Pt"];

// Avalon❗✔️: char *lanthanoids[] =                /* pseudosymbol "lan" */
// Avalon❗✔️:    {
// Avalon❗✔️:        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
// Avalon❗✔️:        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const LANTHANOIDS: &[&str] = &[
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
];

// Avalon❗✔️: char *amino_acids[] =                   /* pseudosymbol "Ami" or "ami"*/
// Avalon❗✔️:    {
// Avalon❗✔️:       "Ala", "Arg", "Asn", "Asp", "Cys",
// Avalon❗✔️:       "Gln", "Glu", "Gly", "His", "Ile",
// Avalon❗✔️:       "Leu", "Lys", "Met", "Phe", "Pro",
// Avalon❗✔️:       "Ser", "Thr", "Trp", "Tyr", "Val",
// Avalon❗✔️:       (char *)NULL,
// Avalon❗✔️:    };
const AMINO_ACIDS: &[&str] = &[
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
    "Thr", "Trp", "Tyr", "Val",
];

fn is_in_string_table(symbol: &str, table: &[&str]) -> bool {
    // Avalon❗✔️: static int IsInStringTable(char *symbol, char *table[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Checks if the string symbol is listed in table[] and returns
    // Avalon❗✔️:  * TRUE if it is and FALSE otherwise.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    char **stringp;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (stringp = table;
    // Avalon❗✔️:         !IsNULL(*stringp);
    // Avalon❗✔️:         stringp++)
    // Avalon❗✔️:       if (0 == strcmp(*stringp,symbol)) return (TRUE);
    // Avalon❗✔️:    return (FALSE);
    // Avalon❗✔️: }
    table.contains(&symbol)
}

fn lower_pseudo_table(token: &str) -> Option<&'static [&'static str]> {
    match token {
        "alk" => Some(ALKALI_METALS),
        "gr2" => Some(GR2),
        "gr3" => Some(GR3),
        "gr4" => Some(GR4),
        "ons" => Some(ONS_TABLE),
        "on2" => Some(ON2),
        "hal" => Some(HALOGENES),
        "ha2" => Some(HA2),
        "trn" => Some(TRANSITION_METALS),
        "tra" => Some(TRA),
        "trb" => Some(TRB),
        "tm1" => Some(TM1),
        "tm2" => Some(TM2),
        "tm3" => Some(TM3),
        "tm4" => Some(TM4),
        "tm5" => Some(TM5),
        "tm6" => Some(TM6),
        "tm7" => Some(TM7),
        "tm8" => Some(TM8),
        "lan" => Some(LANTHANOIDS),
        "ami" => Some(AMINO_ACIDS),
        _ => None,
    }
}

pub(super) fn atom_symbol_match(atom_symbol: &str, pattern: &str) -> bool {
    // Avalon❗✔️: int AtomSymbolMatch(char *atsym, char *pattern)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Returns TRUE if atsym is in the comma delimited list of atom symbols
    // Avalon❗✔️:  * stored in pattern and FALSE otherwise.
    // Avalon❗✔️:  * There are also a number of standard atom type lists like "alk" for alkali metals or
    // Avalon❗✔️:  * "Q" for non-C/non-H defined above as arrays of strings.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    /* static */ char pat_buf[400];
    // Avalon❗✔️:    char *tokp;
    // Avalon❗✔️:
    // Avalon❗✔️:    strcpy(pat_buf,pattern);
    // Avalon❗✔️:    for (tokp = strtok(pat_buf,",");
    // Avalon❗✔️:         !IsNULL(tokp);
    // Avalon❗✔️:         tokp = strtok((char *)NULL,","))
    // Avalon❗✔️:    {
    for token in pattern.split(',').filter(|token| !token.is_empty()) {
        // Avalon❗✔️:       if (islower(*tokp))
        // Avalon❗✔️:       {
        // Avalon❗✔️:          if (0 == strcmp("alk",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,alkali_metals)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("gr2",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,gr2)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("gr3",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,gr3)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("gr4",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,gr4)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("ons",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,ONS_table)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("on2",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,on2)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("hal",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,halogenes)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("ha2",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,ha2)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("trn",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,transition_metals)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tra",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tra)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("trb",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,trb)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm1",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm1)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm2",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm2)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm3",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm3)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm4",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm4)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm5",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm5)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm6",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm6)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm7",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm7)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("tm8",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,tm8)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("lan",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,lanthanoids)) return (TRUE);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else if (0 == strcmp("ami",tokp))
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (IsInStringTable(atsym,amino_acids)) return (TRUE);
        // Avalon❗✔️:          }
        if token.as_bytes().first().is_some_and(u8::is_ascii_lowercase)
            && lower_pseudo_table(token).is_some_and(|table| is_in_string_table(atom_symbol, table))
        {
            return true;
        }
        // Avalon❗✔️:       }
        // Avalon❗✔️:       if (0 == strcmp(atsym,tokp)) return (TRUE);
        if atom_symbol == token {
            return true;
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    if (0 == strcmp("A",pattern))
    // Avalon❗✔️:       return (0 != strcmp("H",atsym));
    if pattern == "A" {
        return atom_symbol != "H";
    }
    // Avalon❗✔️:    else if (0 == strcmp("Qs",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,non_metal_small_solution));
    // Avalon❗✔️:    else if (0 == strcmp("G",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,HC_table));
    // Avalon❗✔️:    else if (0 == strcmp("ONS",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,ONS_table));
    // Avalon❗✔️:    else if (0 == strcmp("X",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,halogenes));
    // Avalon❗✔️:    else if (0 == strcmp("Q",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,non_metal_hetero_elements));
    // Avalon❗✔️:    else if (0 == strcmp("M",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,metals));
    // Avalon❗✔️:    else if (0 == strcmp("Ami",pattern))
    // Avalon❗✔️:       return (IsInStringTable(atsym,amino_acids));
    let whole_pattern_table = match pattern {
        "Qs" => Some(NON_METAL_SMALL_SOLUTION),
        "G" => Some(HC_TABLE),
        "ONS" => Some(ONS_TABLE),
        "X" => Some(HALOGENES),
        "Q" => Some(NON_METAL_HETERO_ELEMENTS),
        "M" => Some(METALS),
        "Ami" => Some(AMINO_ACIDS),
        _ => None,
    };
    if whole_pattern_table.is_some_and(|table| is_in_string_table(atom_symbol, table)) {
        return true;
    }
    // Avalon❗✔️:    return (FALSE);
    // Avalon❗✔️: }
    false
}

pub(super) fn get_symbol_list(molecule: &MoleculeState, index: usize) -> Option<&SymbolList> {
    // Avalon❗✔️: struct symbol_list_t *getList(struct reaccs_molecule_t *mp, int index)
    // Avalon❗✔️: {
    // Avalon❗✔️:    struct symbol_list_t *slp;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (slp = mp->symbol_lists; !IsNULL(slp); slp=slp->next)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       if (slp->atom == index+1) break;
    // Avalon❗✔️:    }
    // Avalon❗✔️:    return (slp);
    // Avalon❗✔️: }
    let atom_number = index.checked_add(1).and_then(|value| i32::try_from(value).ok())?;
    molecule
        .symbol_lists
        .iter()
        .find(|symbol_list| symbol_list.atom == atom_number)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_lowercase_pseudo_symbol_table_has_native_membership_boundaries() {
        let cases = [
            ("alk", "Li", "Be"),
            ("gr2", "Mg", "Al"),
            ("gr3", "Ga", "Ge"),
            ("gr4", "Sn", "Sb"),
            ("ons", "N", "P"),
            ("on2", "Te", "Cl"),
            ("hal", "F", "H"),
            ("ha2", "Cl", "F"),
            ("trn", "Sc", "Al"),
            ("tra", "Pd", "Cu"),
            ("trb", "Ag", "Pd"),
            ("tm1", "Au", "Hg"),
            ("tm2", "Cd", "Ag"),
            ("tm3", "La", "Ce"),
            ("tm4", "Zr", "Nb"),
            ("tm5", "Ta", "W"),
            ("tm6", "Mo", "Tc"),
            ("tm7", "Re", "Os"),
            ("tm8", "Rh", "Ag"),
            ("lan", "Lu", "La"),
            ("ami", "Trp", "C"),
        ];
        for (pattern, member, nonmember) in cases {
            assert!(atom_symbol_match(member, pattern), "{member} in {pattern}");
            assert!(!atom_symbol_match(nonmember, pattern), "{nonmember} not in {pattern}");
        }
    }

    #[test]
    fn whole_pattern_fallbacks_run_after_literal_token_matching() {
        assert!(atom_symbol_match("Xe", "A"));
        assert!(!atom_symbol_match("H", "A"));
        assert!(atom_symbol_match("C", "Qs"));
        assert!(atom_symbol_match("H", "G"));
        assert!(atom_symbol_match("S", "ONS"));
        assert!(atom_symbol_match("Br", "X"));
        assert!(atom_symbol_match("As", "Q"));
        assert!(!atom_symbol_match("C", "Q"));
        assert!(atom_symbol_match("Pu", "M"));
        assert!(!atom_symbol_match("Rn", "M"));
        assert!(atom_symbol_match("Gly", "Ami"));
        assert!(atom_symbol_match("Q", "Q"));
    }

    #[test]
    fn comma_tokenization_matches_strtok_empty_field_and_literal_rules() {
        assert!(atom_symbol_match("N", ",C,,N,"));
        assert!(atom_symbol_match("Li", "C,alk"));
        assert!(atom_symbol_match("alk", "alk"));
        assert!(atom_symbol_match("*", "C,N,A,*"));
        assert!(!atom_symbol_match("O", "C,N,A,*"));
        assert!(!atom_symbol_match("H", ""));
    }

    #[test]
    fn symbol_list_lookup_uses_source_linked_list_order_and_one_based_atoms() {
        let molecule = MoleculeState {
            symbol_lists: vec![
                SymbolList {
                    atom: 2,
                    inclusive: true,
                    symbols: "N,O".to_string(),
                },
                SymbolList {
                    atom: 2,
                    inclusive: false,
                    symbols: "C".to_string(),
                },
            ],
            ..MoleculeState::default()
        };

        let symbol_list = get_symbol_list(&molecule, 1).unwrap();
        assert!(symbol_list.inclusive);
        assert_eq!(symbol_list.symbols, "N,O");
        assert!(get_symbol_list(&molecule, 0).is_none());
    }
}
