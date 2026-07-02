//! RDKit molecular descriptor parity surfaces.
//!
//! The functions in this module intentionally fail closed until their RDKit
//! source paths are ported. They exist so descriptor parity tests can be added
//! before implementation without silently using approximations.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, VecDeque},
    sync::OnceLock,
};

use crate::{
    AdjacencyList, AtomId, BondOrder, Molecule, SubstructMatchParams, SubstructMatchResult,
    chemistry::valence::{
        ValenceModel, assign_valence_with_options, rdkit_atomic_mass, rdkit_element_symbol,
        rdkit_most_common_isotope_mass,
    },
    symmetrize_sssr,
};

const RDKIT_ELECTRON_MASS: f64 = 0.00054857991;
const RDKIT_NUM_HBD_SMARTS: &str = "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]";
const RDKIT_NUM_HBA_SMARTS: &str = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]";
const RDKIT_ROTATABLE_BONDS_NON_STRICT_SMARTS: &str = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
const RDKIT_ROTATABLE_BONDS_STRICT_SMARTS: &str = concat!(
    "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(",
    "[CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=",
    "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$",
    "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([",
    "CH3])[CH3])&!$([CH3])]"
);
const RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_BASE_SMARTS: &str =
    "[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]";
const RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_NON_RING_AMIDES_SMARTS: &str = "[C&!R](=O)NC";
const RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_SYM_RINGS_SMARTS: &str = concat!(
    "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])(",
    "a[!#1])a)]"
);
const RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_TERMINAL_TRIPLE_BONDS_SMARTS: &str = "C#[#6,#7]";
const RDKIT_CRIPPEN_DEFAULT_PARAM_DATA: &str = r#"#ID	SMARTS	logP	MR	Notes/Questions
C1	[CH4]	0.1441	2.503	
C1	[CH3]C	0.1441	2.503	
C1	[CH2](C)C	0.1441	2.503	
C2	[CH](C)(C)C	0	2.433	
C2	[C](C)(C)(C)C	0	2.433	
C3	[CH3][N,O,P,S,F,Cl,Br,I]	-0.2035	2.753	
C3	[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]	-0.2035	2.753	
C4	[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]	-0.2051	2.731	
C4	[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]	-0.2051	2.731	
C5	[C]=[!C;A;!#1]	-0.2783	5.007	
C6	[CH2]=C	0.1551	3.513	
C6	[CH1](=C)[A;!#1]	0.1551	3.513	
C6	[CH0](=C)([A;!#1])[A;!#1]	0.1551	3.513	
C6	[C](=C)=C	0.1551	3.513	
C7	[CX2]#[A;!#1]	0.0017	3.888	
C8	[CH3]c	0.08452	2.464	
C9	[CH3]a	-0.1444	2.412	
C10	[CH2X4]a	-0.0516	2.488	
C11	[CHX4]a	0.1193	2.582	
C12	[CH0X4]a	-0.0967	2.576	
C13	[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]	-0.5443	4.041	
C14	[c][#9]	0	3.257	
C15	[c][#17]	0.245	3.564	
C16	[c][#35]	0.198	3.18	
C17	[c][#53]	0	3.104	
C18	[cH]	0.1581	3.35	
C19	[c](:a)(:a):a	0.2955	4.346	
C20	[c](:a)(:a)-a	0.2713	3.904	
C21	[c](:a)(:a)-C	0.136	3.509	
C22	[c](:a)(:a)-N	0.4619	4.067	
C23	[c](:a)(:a)-O	0.5437	3.853	
C24	[c](:a)(:a)-S	0.1893	2.673	
C25	[c](:a)(:a)=[C,N,O]	-0.8186	3.135	
C26	[C](=C)(a)[A;!#1]	0.264	4.305	
C26	[C](=C)(c)a	0.264	4.305	
C26	[CH1](=C)a	0.264	4.305	
C26	[C]=c	0.264	4.305	
C27	[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]	0.2148	2.693	
CS	[#6]	0.08129	3.243	
H1	[#1][#6,#1]	0.123	1.057	
H2	[#1]O[CX4,c]	-0.2677	1.395	
H2	[#1]O[!#6;!#7;!#8;!#16]	-0.2677	1.395	
H2	[#1][!#6;!#7;!#8]	-0.2677	1.395	
H3	[#1][#7]	0.2142	0.9627	
H3	[#1]O[#7]	0.2142	0.9627	
H4	[#1]OC=[#6,#7,O,S]	0.298	1.805	
H4	[#1]O[O,S]	0.298	1.805	
HS	[#1]	0.1125	1.112	
N1	[NH2+0][A;!#1]	-1.019	2.262	
N2	[NH+0]([A;!#1])[A;!#1]	-0.7096	2.173	
N3	[NH2+0]a	-1.027	2.827	
N4	[NH1+0]([!#1;A,a])a	-0.5188	3	
N5	[NH+0]=[!#1;A,a]	0.08387	1.757	
N6	[N+0](=[!#1;A,a])[!#1;A,a]	0.1836	2.428	
N7	[N+0]([A;!#1])([A;!#1])[A;!#1]	-0.3187	1.839	
N8	[N+0](a)([!#1;A,a])[A;!#1]	-0.4458	2.819	
N8	[N+0](a)(a)a	-0.4458	2.819	
N9	[N+0]#[A;!#1]	0.01508	1.725	
N10	[NH3,NH2,NH;+,+2,+3]	-1.95		
N11	[n+0]	-0.3239	2.202	
N12	[n;+,+2,+3]	-1.119		
N13	[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]	-0.3396	0.2604	
N13	[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]	-0.3396	0.2604	
N13	[NH0;+,+2,+3](=[#6])=[#7]	-0.3396	0.2604	
N14	[N;+,+2,+3]#[A;!#1]	0.2887	3.359	
N14	[N;-,-2,-3]	0.2887	3.359	
N14	[N;+,+2,+3](=[N;-,-2,-3])=N	0.2887	3.359	
NS	[#7]	-0.4806	2.134	
O1	[o]	0.1552	1.08	
O2	[OH,OH2]	-0.2893	0.8238	
O3	[O]([A;!#1])[A;!#1]	-0.0684	1.085	
O4	[O](a)[!#1;A,a]	-0.4195	1.182	
O5	[O]=[#7,#8]	0.0335	3.367	
O5	[OX1;-,-2,-3][#7]	0.0335	3.367	
O6	[OX1;-,-2,-2][#16]	-0.3339	0.7774	
O6	[O;-0]=[#16;-0]	-0.3339	0.7774	
O12	[O-]C(=O)	-1.326		"order flip here intentional"
O7	[OX1;-,-2,-3][!#1;!N;!S]	-1.189	0	
O8	[O]=c	0.1788	3.135	
O9	[O]=[CH]C	-0.1526	0	
O9	[O]=C(C)([A;!#1])	-0.1526	0	
O9	[O]=[CH][N,O]	-0.1526	0	
O9	[O]=[CH2]	-0.1526	0	
O9	[O]=[CX2]=O	-0.1526	0	
O10	[O]=[CH]c	0.1129	0.2215	
O10	[O]=C([C,c])[a;!#1]	0.1129	0.2215	
O10	[O]=C(c)[A;!#1]	0.1129	0.2215	
O11	[O]=C([!#1;!#6])[!#1;!#6]	0.4833	0.389	
OS	[#8]	-0.1188	0.6865	
F	[#9-0]	0.4202	1.108	
Cl	[#17-0]	0.6895	5.853	
Br	[#35-0]	0.8456	8.927	
I	[#53-0]	0.8857	14.02	
Hal	[#9,#17,#35,#53;-]	-2.996		
Hal	[#53;+,+2,+3]	-2.996		
Hal	[+;#3,#11,#19,#37,#55]	-2.996		"Footnote h indicates these should be here?"
P	[#15]	0.8612	6.92	
S2	[S;-,-2,-3,-4,+1,+2,+3,+5,+6]	-0.0024	7.365	"Order flip here is intentional"
S2	[S-0]=[N,O,P,S]	-0.0024	7.365	"Expanded definition of (pseudo-)ionic S"
S1	[S;A]	0.6482	7.591	"Order flip here is intentional"
S3	[s;a]	0.6237	6.691	
Me1	[#3,#11,#19,#37,#55]	-0.3808	5.754	
Me1	[#4,#12,#20,#38,#56]	-0.3808	5.754	
Me1	[#5,#13,#31,#49,#81]	-0.3808	5.754	
Me1	[#14,#32,#50,#82]	-0.3808	5.754	
Me1	[#33,#51,#83]	-0.3808	5.754	
Me1	[#34,#52,#84]	-0.3808	5.754	
Me2	[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]	-0.0025		
Me2	[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]	-0.0025		
Me2	[#72,#73,#74,#75,#76,#77,#78,#79,#80]	-0.0025		
"#;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DescriptorError {
    #[error(
        "RDKit descriptor function {function} is unsupported until {rdkit_function} is source-ported"
    )]
    Unsupported {
        function: &'static str,
        rdkit_function: &'static str,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrippenDescriptorValues {
    pub logp: f64,
    pub molar_refractivity: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumRotatableBondsOptions {
    Default,
    NonStrict,
    Strict,
    StrictLinkages,
}

pub type DescriptorResult<T> = Result<T, DescriptorError>;

pub fn calc_mol_wt(mol: &Molecule, only_heavy: bool) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcAMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getAvgMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: double getAvgMolWt(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   const PeriodicTable *table = PeriodicTable::getTable();
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     if (!onlyHeavy || atom->getAtomicNum() != 1) {
    // RDKit✔️✔️:       res += atom->getMass();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // add our implicit Hs if we need to:
    // RDKit✔️✔️:     if (!onlyHeavy) {
    // RDKit✔️✔️:       res += atom->getTotalNumHs() * table->getAtomicWeight(1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let mut res = 0.0;
    let hydrogen_atomic_weight = rdkit_atomic_mass(1, None);
    let valence =
        assign_valence_with_options(mol, ValenceModel::RdkitLike, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_mol_wt",
                rdkit_function: "Descriptors::calcAMW/MolOps::getAvgMolWt",
            }
        })?;
    for (idx, atom) in mol.atoms().iter().enumerate() {
        if !only_heavy || atom.atomic_number() != 1 {
            res += rdkit_atomic_mass(atom.atomic_number(), atom.isotope());
        }
        if !only_heavy {
            let total_hs = u32::from(atom.explicit_hydrogens())
                + valence
                    .implicit_hydrogens
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
                    .max(0) as u32;
            res += f64::from(total_hs) * hydrogen_atomic_weight;
        }
    }
    Ok(res)
}

pub fn calc_exact_mol_wt(mol: &Molecule, only_heavy: bool) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcExactMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getExactMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: const double electronMass =
    // RDKit✔️✔️:     0.00054857991;  // value from http://physics.nist.gov/cgi-bin/cuu/Value?me
    //
    // RDKit✔️✔️: double getExactMolWt(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   int nHsToCount = 0;
    // RDKit✔️✔️:   const PeriodicTable *table = PeriodicTable::getTable();
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     auto atNum = atom->getAtomicNum();
    // RDKit✔️✔️:     if (atNum != 1 || !onlyHeavy) {
    // RDKit✔️✔️:       if (!atom->getIsotope()) {
    // RDKit✔️✔️:         res += table->getMostCommonIsotopeMass(atNum);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res += atom->getMass();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res -= constants::electronMass * atom->getFormalCharge();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // add our implicit Hs if we need to:
    // RDKit✔️✔️:     if (!onlyHeavy) {
    // RDKit✔️✔️:       nHsToCount += atom->getTotalNumHs(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!onlyHeavy) {
    // RDKit✔️✔️:     res += nHsToCount * table->getMostCommonIsotopeMass(1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let mut res = 0.0;
    let mut hydrogens_to_count = 0_i32;
    let valence =
        assign_valence_with_options(mol, ValenceModel::RdkitLike, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_exact_mol_wt",
                rdkit_function: "Descriptors::calcExactMW/MolOps::getExactMolWt",
            }
        })?;
    for (idx, atom) in mol.atoms().iter().enumerate() {
        let atomic_number = atom.atomic_number();
        if atomic_number != 1 || !only_heavy {
            if atom.isotope().is_none() {
                res += rdkit_most_common_isotope_mass(atomic_number).map_err(|_| {
                    DescriptorError::Unsupported {
                        function: "calc_exact_mol_wt",
                        rdkit_function: "PeriodicTable::getMostCommonIsotopeMass",
                    }
                })?;
            } else {
                res += rdkit_atomic_mass(atomic_number, atom.isotope());
            }
            res -= RDKIT_ELECTRON_MASS * f64::from(atom.formal_charge());
        }
        if !only_heavy {
            hydrogens_to_count += i32::from(atom.explicit_hydrogens())
                + valence
                    .implicit_hydrogens
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
                    .max(0);
        }
    }
    if !only_heavy {
        res += f64::from(hydrogens_to_count)
            * rdkit_most_common_isotope_mass(1).map_err(|_| DescriptorError::Unsupported {
                function: "calc_exact_mol_wt",
                rdkit_function: "PeriodicTable::getMostCommonIsotopeMass",
            })?;
    }
    Ok(res)
}

pub fn calc_mol_formula(
    mol: &Molecule,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
) -> DescriptorResult<String> {
    // RDKit✔️✔️: std::string calcMolFormula(const ROMol &mol, bool separateIsotopes,
    // RDKit✔️✔️:                              bool abbreviateHIsotopes) {
    // RDKit✔️✔️:   return MolOps::getMolFormula(mol, separateIsotopes, abbreviateHIsotopes);
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: std::string getMolFormula(const ROMol &mol, bool separateIsotopes,
    // RDKit✔️✔️:                           bool abbreviateHIsotopes) {
    // RDKit✔️✔️:   std::ostringstream res;
    // RDKit✔️✔️:   std::map<std::pair<unsigned int, std::string>, unsigned int> counts;
    // RDKit✔️✔️:   int charge = 0;
    // RDKit✔️✔️:   unsigned int nHs = 0;
    // RDKit✔️✔️:   const PeriodicTable *table = PeriodicTable::getTable();
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     int atNum = atom->getAtomicNum();
    // RDKit✔️✔️:     std::pair<unsigned int, std::string> key =
    // RDKit✔️✔️:         std::make_pair(0, table->getElementSymbol(atNum));
    // RDKit✔️✔️:     if (separateIsotopes) {
    // RDKit✔️✔️:       unsigned int isotope = atom->getIsotope();
    // RDKit✔️✔️:       if (abbreviateHIsotopes && atNum == 1 && (isotope == 2 || isotope == 3)) {
    // RDKit✔️✔️:         if (isotope == 2) {
    // RDKit✔️✔️:           key.second = "D";
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           key.second = "T";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         key.first = isotope;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     counts[key] += 1;
    // RDKit✔️✔️:     nHs += atom->getTotalNumHs();
    // RDKit✔️✔️:     charge += atom->getFormalCharge();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (nHs) {
    // RDKit✔️✔️:     counts[{0, "H"}] += nHs;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::list<std::pair<unsigned int, std::string>> ks;
    // RDKit✔️✔️:   for (const auto &pr : counts) {
    // RDKit✔️✔️:     ks.push_back(pr.first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ks.sort(HillCompare);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &key : ks) {
    // RDKit✔️✔️:     if (key.first > 0) {
    // RDKit✔️✔️:       res << "[" << key.first << key.second << "]";
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res << key.second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (counts[key] > 1) {
    // RDKit✔️✔️:       res << counts[key];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (charge > 0) {
    // RDKit✔️✔️:     res << "+";
    // RDKit✔️✔️:     if (charge > 1) {
    // RDKit✔️✔️:       res << charge;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (charge < 0) {
    // RDKit✔️✔️:     res << "-";
    // RDKit✔️✔️:     if (charge < -1) {
    // RDKit✔️✔️:       res << -1 * charge;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res.str();
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let valence =
        assign_valence_with_options(mol, ValenceModel::RdkitLike, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_mol_formula",
                rdkit_function: "Atom::getTotalNumHs",
            }
        })?;
    let mut counts = BTreeMap::<FormulaKey, u32>::new();
    let mut charge = 0_i32;
    let mut hydrogens = 0_u32;
    for (idx, atom) in mol.atoms().iter().enumerate() {
        let atomic_number = atom.atomic_number();
        let mut key = FormulaKey {
            isotope: 0,
            symbol: rdkit_element_symbol(atomic_number).map_err(|_| {
                DescriptorError::Unsupported {
                    function: "calc_mol_formula",
                    rdkit_function: "PeriodicTable::getElementSymbol",
                }
            })?,
        };
        if separate_isotopes {
            let isotope = atom.isotope().map(u32::from).unwrap_or(0);
            if abbreviate_h_isotopes && atomic_number == 1 && (isotope == 2 || isotope == 3) {
                if isotope == 2 {
                    key.symbol = "D";
                } else {
                    key.symbol = "T";
                }
            } else {
                key.isotope = isotope;
            }
        }
        *counts.entry(key).or_insert(0) += 1;
        hydrogens += rdkit_total_num_hs_without_neighbors(
            atom.explicit_hydrogens(),
            valence.implicit_hydrogens.get(idx).copied().unwrap_or(0),
        )?;
        charge += i32::from(atom.formal_charge());
    }
    if hydrogens != 0 {
        *counts
            .entry(FormulaKey {
                isotope: 0,
                symbol: "H",
            })
            .or_insert(0) += hydrogens;
    }
    let mut keys = counts.keys().copied().collect::<Vec<_>>();
    keys.sort_by(hill_compare_ordering);

    let mut res = String::new();
    for key in keys {
        if key.isotope > 0 {
            res.push('[');
            res.push_str(&key.isotope.to_string());
            res.push_str(key.symbol);
            res.push(']');
        } else {
            res.push_str(key.symbol);
        }
        let count = counts[&key];
        if count > 1 {
            res.push_str(&count.to_string());
        }
    }
    if charge > 0 {
        res.push('+');
        if charge > 1 {
            res.push_str(&charge.to_string());
        }
    } else if charge < 0 {
        res.push('-');
        if charge < -1 {
            res.push_str(&(-charge).to_string());
        }
    }
    Ok(res)
}

pub fn calc_num_hbd(mol: &Molecule) -> DescriptorResult<u32> {
    // RDKit✔️✔️: class ss_matcher {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    // RDKit✔️✔️:     m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    // RDKit✔️✔️:     RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit✔️✔️:     m_matcher = p;
    // RDKit✔️✔️:     POSTCONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   const RDKit::ROMol *getMatcher() const { return m_matcher; };
    // RDKit✔️✔️:   unsigned int countMatches(const RDKit::ROMol &mol) const {
    // RDKit✔️✔️:     PRECONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:     std::vector<RDKit::MatchVectType> matches;
    // RDKit✔️✔️:     // This is an ugly one. Recursive queries aren't thread safe.
    // RDKit✔️✔️:     // Unfortunately we have to take a performance hit here in order
    // RDKit✔️✔️:     // to guarantee thread safety
    // RDKit✔️✔️:     if (m_needCopies) {
    // RDKit✔️✔️:       const RDKit::ROMol nm(*(m_matcher), true);
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       const RDKit::ROMol &nm = *m_matcher;
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return matches.size();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ~ss_matcher() { delete m_matcher; };
    // RDKit✔️✔️:
    // RDKit✔️✔️:  private:
    // RDKit✔️✔️:   ss_matcher() : m_pattern("") {};
    // RDKit✔️✔️:   std::string m_pattern;
    // RDKit✔️✔️:   bool m_needCopies{false};
    // RDKit✔️✔️:   const RDKit::ROMol *m_matcher{nullptr};
    // RDKit✔️✔️: };
    //
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    //
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHBD, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]",
    // RDKit✔️✔️:                 "2.0.1");
    rdkit_count_smarts_matches("calc_num_hbd", RDKIT_NUM_HBD_SMARTS)
        .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_hbd"))
}

pub fn calc_num_hba(mol: &Molecule) -> DescriptorResult<u32> {
    // RDKit✔️✔️: class ss_matcher {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    // RDKit✔️✔️:     m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    // RDKit✔️✔️:     RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit✔️✔️:     m_matcher = p;
    // RDKit✔️✔️:     POSTCONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   const RDKit::ROMol *getMatcher() const { return m_matcher; };
    // RDKit✔️✔️:   unsigned int countMatches(const RDKit::ROMol &mol) const {
    // RDKit✔️✔️:     PRECONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:     std::vector<RDKit::MatchVectType> matches;
    // RDKit✔️✔️:     // This is an ugly one. Recursive queries aren't thread safe.
    // RDKit✔️✔️:     // Unfortunately we have to take a performance hit here in order
    // RDKit✔️✔️:     // to guarantee thread safety
    // RDKit✔️✔️:     if (m_needCopies) {
    // RDKit✔️✔️:       const RDKit::ROMol nm(*(m_matcher), true);
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       const RDKit::ROMol &nm = *m_matcher;
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return matches.size();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ~ss_matcher() { delete m_matcher; };
    // RDKit✔️✔️:
    // RDKit✔️✔️:  private:
    // RDKit✔️✔️:   ss_matcher() : m_pattern("") {};
    // RDKit✔️✔️:   std::string m_pattern;
    // RDKit✔️✔️:   bool m_needCopies{false};
    // RDKit✔️✔️:   const RDKit::ROMol *m_matcher{nullptr};
    // RDKit✔️✔️: };
    //
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    //
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHBA,
    // RDKit✔️✔️:                 "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$("
    // RDKit✔️✔️:                 "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]",
    // RDKit✔️✔️:                 "2.0.2");
    rdkit_count_smarts_matches("calc_num_hba", RDKIT_NUM_HBA_SMARTS)
        .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_hba"))
}

pub fn calc_fraction_csp3(mol: &Molecule) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcFractionCSP3(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int nCSP3 = 0;
    // RDKit✔️✔️:   unsigned int nC = 0;
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
    // RDKit✔️✔️:     if (at->getAtomicNum() == 6) {
    // RDKit✔️✔️:       ++nC;
    // RDKit✔️✔️:       if (at->getTotalDegree() == 4) {
    // RDKit✔️✔️:         ++nCSP3;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atBegin;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!nC) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<double>(nCSP3) / nC;
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: unsigned int Atom::getDegree() const {
    // RDKit✔️✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let valence =
        assign_valence_with_options(mol, ValenceModel::RdkitLike, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_fraction_csp3",
                rdkit_function: "Atom::getTotalDegree",
            }
        })?;
    let adjacency = &mol.topology_block().adjacency;
    let mut n_csp3 = 0_u32;
    let mut n_c = 0_u32;
    for (idx, atom) in mol.atoms().iter().enumerate() {
        if atom.atomic_number() == 6 {
            n_c += 1;
            let attached_hs = u32::from(atom.explicit_hydrogens())
                + valence
                    .implicit_hydrogens
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
                    .max(0) as u32;
            let total_degree = adjacency.neighbors_of(idx).len() as u32 + attached_hs;
            if total_degree == 4 {
                n_csp3 += 1;
            }
        }
    }
    if n_c == 0 {
        return Ok(0.0);
    }
    Ok(f64::from(n_csp3) / f64::from(n_c))
}

pub fn calc_crippen_descriptors(
    mol: &Molecule,
    include_hs: bool,
    _force: bool,
) -> DescriptorResult<CrippenDescriptorValues> {
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_crippenLogP)) {
    // RDKit✔️✔️:     mol.getProp(common_properties::_crippenLogP, logp);
    // RDKit✔️✔️:     mol.getProp(common_properties::_crippenMR, mr);
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // COSMolKit modeled input state: RDKit's mutable typed computed-property
    // cache entries `_crippenLogP` and `_crippenMR` are not modeled by this
    // value-style `&Molecule` descriptor API, so no cached descriptor input can
    // be observed or reused here.
    //
    // RDKit✔️❌:   // this isn't as bad as it looks, we aren't actually going
    // RDKit✔️❌:   // to harm the molecule in any way!
    // RDKit✔️❌:   auto *workMol = const_cast<ROMol *>(&mol);
    // RDKit✔️❌:   if (includeHs) {
    // RDKit✔️❌:     workMol = MolOps::addHs(mol, false, false);
    // RDKit✔️❌:   }
    let work_mol;
    let work_ref = if include_hs {
        work_mol = mol
            .with_hydrogens()
            .map_err(|_| DescriptorError::Unsupported {
                function: "calc_crippen_descriptors",
                rdkit_function: "MolOps::addHs(mol, false, false)",
            })?;
        &work_mol
    } else {
        mol
    };
    // RDKit✔️✔️:   std::vector<double> logpContribs(workMol->getNumAtoms());
    // RDKit✔️✔️:   std::vector<double> mrContribs(workMol->getNumAtoms());
    // RDKit✔️✔️:   getCrippenAtomContribs(*workMol, logpContribs, mrContribs, force);
    let contribs = rdkit_crippen_atom_contribs(work_ref)?;
    // RDKit✔️✔️:   logp = 0.0;
    // RDKit✔️✔️:   for (std::vector<double>::const_iterator iter = logpContribs.begin();
    // RDKit✔️✔️:        iter != logpContribs.end(); ++iter) {
    // RDKit✔️✔️:     logp += *iter;
    // RDKit✔️✔️:   }
    let logp = contribs.logp.iter().copied().sum();
    // RDKit✔️✔️:   mr = 0.0;
    // RDKit✔️✔️:   for (std::vector<double>::const_iterator iter = mrContribs.begin();
    // RDKit✔️✔️:        iter != mrContribs.end(); ++iter) {
    // RDKit✔️✔️:     mr += *iter;
    // RDKit✔️✔️:   }
    let molar_refractivity = contribs.mr.iter().copied().sum();
    //
    // RDKit✔️✔️:   if (includeHs) {
    // RDKit✔️✔️:     delete workMol;
    // RDKit✔️✔️:   }
    //
    // RDKit✔️✔️:   mol.setProp(common_properties::_crippenLogP, logp, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_crippenMR, mr, true);
    // COSMolKit modeled input state: Crippen descriptor calls return
    // `CrippenDescriptorValues` and do not mutate the input molecule with RDKit
    // computed-property cache entries.
    // RDKit✔️✔️: };
    Ok(CrippenDescriptorValues {
        logp,
        molar_refractivity,
    })
}

#[derive(Debug, Clone)]
struct CrippenParam {
    idx: usize,
    label: &'static str,
    smarts: &'static str,
    logp: f64,
    mr: f64,
}

fn rdkit_crippen_params() -> &'static [CrippenParam] {
    static PARAMS: OnceLock<Vec<CrippenParam>> = OnceLock::new();
    PARAMS.get_or_init(rdkit_parse_default_crippen_params)
}

fn rdkit_parse_default_crippen_params() -> Vec<CrippenParam> {
    // RDKit✔️✔️: CrippenParamCollection::CrippenParamCollection(const std::string &paramData) {
    // RDKit✔️✔️:   std::string params;
    // RDKit✔️✔️:   boost::char_separator<char> tabSep("\t", "", boost::keep_empty_tokens);
    // RDKit✔️✔️:   if (paramData == "") {
    // RDKit✔️✔️:     params = defaultParamData;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     params = paramData;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::istringstream inStream(params);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (!inStream.eof()) {
    // RDKit✔️✔️:     if (inLine[0] != '#') {
    // RDKit✔️✔️:       CrippenParams paramObj;
    // RDKit✔️✔️:       paramObj.idx = idx++;
    // RDKit✔️✔️:       tokenizer tokens(inLine, tabSep);
    // RDKit✔️✔️:       tokenizer::iterator token = tokens.begin();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       paramObj.label = *token;
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       paramObj.smarts = *token;
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       if (*token != "") {
    // RDKit✔️✔️:         paramObj.logp = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         paramObj.logp = 0.0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       if (*token != "") {
    // RDKit✔️✔️:         try {
    // RDKit✔️✔️:           paramObj.mr = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:         } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:           paramObj.mr = 0.0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         paramObj.mr = 0.0;
    // RDKit✔️✔️:       }
    // RDKit✔️❌:       paramObj.dp_pattern =
    // RDKit✔️❌:           boost::shared_ptr<const ROMol>(SmartsToMol(paramObj.smarts));
    // COSMolKit preserves descriptor semantics by storing the source SMARTS and
    // compiling it at match time in `rdkit_crippen_atom_contribs()`. This is
    // materially worse than RDKit's cached parsed query molecule but does not
    // change ordered first-match atom typing for the modeled input state.
    // RDKit✔️✔️:       d_params.push_back(paramObj);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut params = Vec::new();
    for line in RDKIT_CRIPPEN_DEFAULT_PARAM_DATA.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut tokens = line.split('\t');
        let Some(label) = tokens.next() else {
            continue;
        };
        let Some(smarts) = tokens.next() else {
            continue;
        };
        let logp = tokens
            .next()
            .filter(|value| !value.is_empty())
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        let mr = tokens
            .next()
            .filter(|value| !value.is_empty())
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        params.push(CrippenParam {
            idx: params.len(),
            label,
            smarts,
            logp,
            mr,
        });
    }
    params
}

#[derive(Debug, Clone)]
struct CrippenAtomContribs {
    logp: Vec<f64>,
    mr: Vec<f64>,
    atom_types: Vec<usize>,
    atom_type_labels: Vec<&'static str>,
}

fn rdkit_crippen_atom_contribs(mol: &Molecule) -> DescriptorResult<CrippenAtomContribs> {
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_crippenLogPContribs)) {
    // RDKit✔️✔️:     std::vector<double> tmpVect1, tmpVect2;
    // RDKit✔️✔️:     mol.getProp(common_properties::_crippenLogPContribs, tmpVect1);
    // RDKit✔️✔️:     mol.getProp(common_properties::_crippenMRContribs, tmpVect2);
    // RDKit✔️✔️:     if (tmpVect1.size() == mol.getNumAtoms() &&
    // RDKit✔️✔️:         tmpVect2.size() == mol.getNumAtoms()) {
    // RDKit✔️✔️:       logpContribs = tmpVect1;
    // RDKit✔️✔️:       mrContribs = tmpVect2;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // COSMolKit modeled input state: RDKit's vector-valued computed-property
    // cache entries `_crippenLogPContribs` and `_crippenMRContribs` are not
    // modeled, so cached per-atom contribution input cannot be observed here.
    //
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomNeeded(mol.getNumAtoms());
    // RDKit✔️✔️:   atomNeeded.set();
    // RDKit✔️✔️:   const CrippenParamCollection *params = CrippenParamCollection::getParams();
    // RDKit✔️✔️:   for (const auto &param : *params) {
    // RDKit✔️✔️:     std::vector<MatchVectType> matches;
    // RDKit✔️✔️:     SubstructMatch(mol, *(param.dp_pattern.get()), matches, false, true);
    // RDKit✔️✔️:     for (std::vector<MatchVectType>::const_iterator matchIt = matches.begin();
    // RDKit✔️✔️:          matchIt != matches.end(); ++matchIt) {
    // RDKit✔️✔️:       int idx = (*matchIt)[0].second;
    // RDKit✔️✔️:       if (atomNeeded[idx]) {
    // RDKit✔️✔️:         atomNeeded[idx] = 0;
    // RDKit✔️✔️:         logpContribs[idx] = param.logp;
    // RDKit✔️✔️:         mrContribs[idx] = param.mr;
    // RDKit✔️✔️:         if (atomTypes) {
    // RDKit✔️✔️:           (*atomTypes)[idx] = param.idx;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (atomTypeLabels) {
    // RDKit✔️✔️:           (*atomTypeLabels)[idx] = param.label;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atomNeeded.none()) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setProp(common_properties::_crippenLogPContribs, logpContribs, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_crippenMRContribs, mrContribs, true);
    // COSMolKit modeled input state: per-atom Crippen contributions remain a
    // local private helper result and are not written into RDKit-style computed
    // property cache entries.
    let mut contribs = CrippenAtomContribs {
        logp: vec![0.0; mol.num_atoms()],
        mr: vec![0.0; mol.num_atoms()],
        atom_types: vec![usize::MAX; mol.num_atoms()],
        atom_type_labels: vec![""; mol.num_atoms()],
    };
    let mut atom_needed = vec![true; mol.num_atoms()];
    let params = rdkit_crippen_params();
    let match_params = SubstructMatchParams {
        max_matches: 1000,
        uniquify: false,
        use_chirality: false,
        specified_stereo_query_matches_unspecified: false,
    };
    for param in params {
        let query = rdkit_count_smarts_matches("calc_crippen_descriptors", param.smarts)?;
        let matches = crate::get_substruct_matches_with_params(mol, &query, &match_params);
        for matched in matches {
            let Some(&idx) = matched.atom_mapping.first() else {
                continue;
            };
            if idx < atom_needed.len() && atom_needed[idx] {
                atom_needed[idx] = false;
                contribs.logp[idx] = param.logp;
                contribs.mr[idx] = param.mr;
                contribs.atom_types[idx] = param.idx;
                contribs.atom_type_labels[idx] = param.label;
            }
        }
        if atom_needed.iter().all(|needed| !*needed) {
            break;
        }
    }
    Ok(contribs)
}

pub fn calc_tpsa(mol: &Molecule, force: bool, include_sandp: bool) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcTPSA(const ROMol &mol, bool force, bool includeSandP) {
    // RDKit✔️✔️:   std::string pname =
    // RDKit✔️✔️:       (boost::format("%s-%s") % common_properties::_tpsa % includeSandP).str();
    // RDKit✔️✔️:   if (!force && mol.hasProp(pname)) {
    // RDKit✔️✔️:     double res;
    // RDKit✔️✔️:     mol.getProp(pname, res);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // COSMolKit modeled input state: RDKit's typed computed-property cache
    // entries `_tpsa-true` and `_tpsa-false` are not modeled by this
    // value-style `&Molecule` descriptor API, so cached TPSA input cannot be
    // observed or reused here.
    // RDKit✔️✔️:   std::vector<double> contribs;
    // RDKit✔️✔️:   contribs.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   res = getTPSAAtomContribs(mol, contribs, force, includeSandP);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(rdkit_tpsa_atom_contribs(mol, force, include_sandp)?.total)
}

#[derive(Debug, Clone)]
struct TpsaAtomContribs {
    total: f64,
    atom_contribs: Vec<f64>,
}

fn rdkit_tpsa_atom_contribs(
    mol: &Molecule,
    _force: bool,
    include_sandp: bool,
) -> DescriptorResult<TpsaAtomContribs> {
    // RDKit✔️✔️: double getTPSAAtomContribs(const ROMol &mol, std::vector<double> &Vi,
    // RDKit✔️✔️:                            bool force, bool includeSandP) {
    // RDKit✔️✔️:   TEST_ASSERT(Vi.size() >= mol.getNumAtoms());
    // RDKit✔️✔️:   double res = 0;
    // RDKit✔️✔️:   std::string pname =
    // RDKit✔️✔️:       (boost::format("%s-%s") % common_properties::_tpsa % includeSandP).str();
    // RDKit✔️✔️:   std::string contribsName =
    // RDKit✔️✔️:       (boost::format("%s-%s") % common_properties::_tpsaAtomContribs %
    // RDKit✔️✔️:        includeSandP)
    // RDKit✔️✔️:           .str();
    // RDKit✔️✔️:   if (!force && mol.hasProp(contribsName)) {
    // RDKit✔️✔️:     mol.getProp(contribsName, Vi);
    // RDKit✔️✔️:     mol.getProp(pname, res);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // COSMolKit modeled input state: RDKit's `_tpsaAtomContribs-*` vector
    // cache and paired `_tpsa-*` total cache are not modeled. Per-atom TPSA
    // contributions are private local helper output.
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   std::vector<int> nNbrs(nAtoms, 0), nSing(nAtoms, 0), nDoub(nAtoms, 0),
    // RDKit✔️✔️:       nTrip(nAtoms, 0), nArom(nAtoms, 0), nHs(nAtoms, 0);
    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bIt = mol.beginBonds(); bIt != mol.endBonds();
    // RDKit✔️✔️:        ++bIt) {
    // RDKit✔️✔️:     const Bond *bnd = (*bIt);
    // RDKit✔️✔️:     if (bnd->getBeginAtom()->getAtomicNum() == 1) {
    // RDKit✔️✔️:       nNbrs[bnd->getEndAtomIdx()] -= 1;
    // RDKit✔️✔️:       nHs[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:     } else if (bnd->getEndAtom()->getAtomicNum() == 1) {
    // RDKit✔️✔️:       nNbrs[bnd->getBeginAtomIdx()] -= 1;
    // RDKit✔️✔️:       nHs[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:     } else if (bnd->getIsAromatic()) {
    // RDKit✔️✔️:       nArom[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:       nArom[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       switch (bnd->getBondType()) {
    // RDKit✔️✔️:         case Bond::SINGLE:
    // RDKit✔️✔️:           nSing[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nSing[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::DOUBLE:
    // RDKit✔️✔️:           nDoub[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nDoub[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::TRIPLE:
    // RDKit✔️✔️:           nTrip[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nTrip[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let n_atoms = mol.num_atoms();
    let mut n_nbrs = vec![0_i32; n_atoms];
    let mut n_sing = vec![0_i32; n_atoms];
    let mut n_doub = vec![0_i32; n_atoms];
    let mut n_trip = vec![0_i32; n_atoms];
    let mut n_arom = vec![0_i32; n_atoms];
    let mut n_hs = vec![0_i32; n_atoms];
    for bnd in mol.bonds() {
        let begin_idx = bnd.begin().index();
        let end_idx = bnd.end().index();
        let begin_atom = &mol.atoms()[begin_idx];
        let end_atom = &mol.atoms()[end_idx];
        if begin_atom.atomic_number() == 1 {
            n_nbrs[end_idx] -= 1;
            n_hs[end_idx] += 1;
        } else if end_atom.atomic_number() == 1 {
            n_nbrs[begin_idx] -= 1;
            n_hs[begin_idx] += 1;
        } else if bnd.is_aromatic() {
            n_arom[begin_idx] += 1;
            n_arom[end_idx] += 1;
        } else {
            match bnd.order() {
                BondOrder::Single => {
                    n_sing[begin_idx] += 1;
                    n_sing[end_idx] += 1;
                }
                BondOrder::Double => {
                    n_doub[begin_idx] += 1;
                    n_doub[end_idx] += 1;
                }
                BondOrder::Triple => {
                    n_trip[begin_idx] += 1;
                    n_trip[end_idx] += 1;
                }
                _ => {}
            }
        }
    }

    let rings = symmetrize_sssr(mol).map_err(|_| DescriptorError::Unsupported {
        function: "calc_tpsa",
        rdkit_function: "RingInfo::isAtomInRingOfSize",
    })?;
    let valence =
        assign_valence_with_options(mol, ValenceModel::RdkitLike, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_tpsa",
                rdkit_function: "Atom::getTotalNumHs",
            }
        })?;
    let adjacency = &mol.topology_block().adjacency;
    let mut atom_contribs = vec![0.0; n_atoms];
    let mut res = 0.0;

    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:     int atNum = atom->getAtomicNum();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (atNum != 7 && atNum != 8 &&
    // RDKit✔️✔️:         (!includeSandP || (atNum != 15 && atNum != 16))) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     nHs[i] += atom->getTotalNumHs();
    // RDKit✔️✔️:     int chg = atom->getFormalCharge();
    // RDKit✔️✔️:     bool in3Ring = mol.getRingInfo()->isAtomInRingOfSize(i, 3);
    // RDKit✔️✔️:     nNbrs[i] += atom->getDegree();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     double tmp = -1;
    for i in 0..n_atoms {
        let atom = &mol.atoms()[i];
        let at_num = atom.atomic_number();
        if at_num != 7 && at_num != 8 && (!include_sandp || (at_num != 15 && at_num != 16)) {
            continue;
        }
        n_hs[i] += rdkit_tpsa_total_num_hs_without_neighbors(
            atom.explicit_hydrogens(),
            valence.implicit_hydrogens.get(i).copied().unwrap_or(0),
        )? as i32;
        let chg = i32::from(atom.formal_charge());
        let in3_ring = rings.is_atom_in_ring_of_size(atom.id(), 3);
        n_nbrs[i] += adjacency.neighbors_of(i).len() as i32;
        let mut tmp = -1.0;

        // RDKit✔️✔️:     if (atNum == 7) {
        if at_num == 7 {
            tmp = rdkit_tpsa_nitrogen_contrib(
                n_nbrs[i], n_sing[i], n_doub[i], n_trip[i], n_arom[i], n_hs[i], chg, in3_ring,
            );
        // RDKit✔️✔️:     } else if (atNum == 8) {
        } else if at_num == 8 {
            tmp = rdkit_tpsa_oxygen_contrib(
                n_nbrs[i], n_sing[i], n_doub[i], n_arom[i], n_hs[i], chg, in3_ring,
            );
        // RDKit✔️✔️:     } else if (includeSandP && atNum == 15) {
        } else if include_sandp && at_num == 15 {
            tmp = rdkit_tpsa_phosphorus_contrib(n_nbrs[i], n_sing[i], n_doub[i], n_hs[i], chg);
        // RDKit✔️✔️:     } else if (includeSandP && atNum == 16) {
        } else if include_sandp && at_num == 16 {
            tmp =
                rdkit_tpsa_sulfur_contrib(n_nbrs[i], n_sing[i], n_doub[i], n_arom[i], n_hs[i], chg);
        }
        // RDKit✔️✔️:     Vi[i] = tmp;
        // RDKit✔️✔️:     res += tmp;
        atom_contribs[i] = tmp;
        res += tmp;
    }
    // RDKit✔️✔️:   mol.setProp(contribsName, Vi, true);
    // RDKit✔️✔️:   mol.setProp(pname, res, true);
    // COSMolKit modeled input state: TPSA calls return values and do not write
    // RDKit-style includeSandP-keyed computed-property cache entries into the
    // input molecule.
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(TpsaAtomContribs {
        total: res,
        atom_contribs,
    })
}

#[allow(clippy::too_many_arguments)]
fn rdkit_tpsa_nitrogen_contrib(
    n_nbrs: i32,
    n_sing: i32,
    n_doub: i32,
    n_trip: i32,
    n_arom: i32,
    n_hs: i32,
    chg: i32,
    in3_ring: bool,
) -> f64 {
    // RDKit✔️✔️:     if (atNum == 7) {
    // RDKit✔️✔️:       switch (nNbrs[i]) {
    // RDKit✔️✔️:         case 1:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nTrip[i] == 1) {
    // RDKit✔️✔️:             tmp = 23.79;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 23.85;
    // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 0 && nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 26.02;
    // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 1 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 25.59;
    // RDKit✔️✔️:           } else if (nHs[i] == 3 && chg == 1 && nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 27.64;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 12.36;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nTrip[i] == 1 &&
    // RDKit✔️✔️:                      nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 13.60;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 && in3Ring) {
    // RDKit✔️✔️:             tmp = 21.94;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 && !in3Ring) {
    // RDKit✔️✔️:             tmp = 12.03;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nTrip[i] == 1 &&
    // RDKit✔️✔️:                      nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 4.36;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nDoub[i] == 1 &&
    // RDKit✔️✔️:                      nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 13.97;
    // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 1 && nSing[i] == 2) {
    // RDKit✔️✔️:             tmp = 16.61;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 12.89;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 15.79;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 14.14;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && in3Ring) {
    // RDKit✔️✔️:             tmp = 3.01;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && !in3Ring) {
    // RDKit✔️✔️:             tmp = 3.24;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 &&
    // RDKit✔️✔️:                      nDoub[i] == 2) {
    // RDKit✔️✔️:             tmp = 11.68;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nSing[i] == 2 &&
    // RDKit✔️✔️:                      nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 3.01;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nSing[i] == 3) {
    // RDKit✔️✔️:             tmp = 4.44;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 3) {
    // RDKit✔️✔️:             tmp = 4.41;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 &&
    // RDKit✔️✔️:                      nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 4.93;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1 &&
    // RDKit✔️✔️:                      nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 8.39;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nArom[i] == 3) {
    // RDKit✔️✔️:             tmp = 4.10;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nSing[i] == 1 &&
    // RDKit✔️✔️:                      nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 3.88;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           if (nHs[i] == 0 && nSing[i] == 4 && chg == 1) {
    // RDKit✔️✔️:             tmp = 0.0;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (tmp < 0.0) {
    // RDKit✔️✔️:         tmp = 30.5 - nNbrs[i] * 8.2 + nHs[i] * 1.5;
    // RDKit✔️✔️:         if (tmp < 0) {
    // RDKit✔️✔️:           tmp = 0.0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    let mut tmp = -1.0;
    match n_nbrs {
        1 => {
            if n_hs == 0 && chg == 0 && n_trip == 1 {
                tmp = 23.79;
            } else if n_hs == 1 && chg == 0 && n_doub == 1 {
                tmp = 23.85;
            } else if n_hs == 2 && chg == 0 && n_sing == 1 {
                tmp = 26.02;
            } else if n_hs == 2 && chg == 1 && n_doub == 1 {
                tmp = 25.59;
            } else if n_hs == 3 && chg == 1 && n_sing == 1 {
                tmp = 27.64;
            }
        }
        2 => {
            if n_hs == 0 && chg == 0 && n_sing == 1 && n_doub == 1 {
                tmp = 12.36;
            } else if n_hs == 0 && chg == 0 && n_trip == 1 && n_doub == 1 {
                tmp = 13.60;
            } else if n_hs == 1 && chg == 0 && n_sing == 2 && in3_ring {
                tmp = 21.94;
            } else if n_hs == 1 && chg == 0 && n_sing == 2 && !in3_ring {
                tmp = 12.03;
            } else if n_hs == 0 && chg == 1 && n_trip == 1 && n_sing == 1 {
                tmp = 4.36;
            } else if n_hs == 1 && chg == 1 && n_doub == 1 && n_sing == 1 {
                tmp = 13.97;
            } else if n_hs == 2 && chg == 1 && n_sing == 2 {
                tmp = 16.61;
            } else if n_hs == 0 && chg == 0 && n_arom == 2 {
                tmp = 12.89;
            } else if n_hs == 1 && chg == 0 && n_arom == 2 {
                tmp = 15.79;
            } else if n_hs == 1 && chg == 1 && n_arom == 2 {
                tmp = 14.14;
            }
        }
        3 => {
            if n_hs == 0 && chg == 0 && n_sing == 3 && in3_ring {
                tmp = 3.01;
            } else if n_hs == 0 && chg == 0 && n_sing == 3 && !in3_ring {
                tmp = 3.24;
            } else if n_hs == 0 && chg == 0 && n_sing == 1 && n_doub == 2 {
                tmp = 11.68;
            } else if n_hs == 0 && chg == 1 && n_sing == 2 && n_doub == 1 {
                tmp = 3.01;
            } else if n_hs == 1 && chg == 1 && n_sing == 3 {
                tmp = 4.44;
            } else if n_hs == 0 && chg == 0 && n_arom == 3 {
                tmp = 4.41;
            } else if n_hs == 0 && chg == 0 && n_sing == 1 && n_arom == 2 {
                tmp = 4.93;
            } else if n_hs == 0 && chg == 0 && n_doub == 1 && n_arom == 2 {
                tmp = 8.39;
            } else if n_hs == 0 && chg == 1 && n_arom == 3 {
                tmp = 4.10;
            } else if n_hs == 0 && chg == 1 && n_sing == 1 && n_arom == 2 {
                tmp = 3.88;
            }
        }
        4 => {
            if n_hs == 0 && n_sing == 4 && chg == 1 {
                tmp = 0.0;
            }
        }
        _ => {}
    }
    if tmp < 0.0 {
        tmp = 30.5 - f64::from(n_nbrs) * 8.2 + f64::from(n_hs) * 1.5;
        if tmp < 0.0 {
            tmp = 0.0;
        }
    }
    tmp
}

fn rdkit_tpsa_oxygen_contrib(
    n_nbrs: i32,
    n_sing: i32,
    n_doub: i32,
    n_arom: i32,
    n_hs: i32,
    chg: i32,
    in3_ring: bool,
) -> f64 {
    // RDKit✔️✔️:     } else if (atNum == 8) {
    // RDKit✔️✔️:       switch (nNbrs[i]) {
    // RDKit✔️✔️:         case 1:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 17.07;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 20.23;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == -1 && nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 23.06;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && in3Ring) {
    // RDKit✔️✔️:             tmp = 12.53;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && !in3Ring) {
    // RDKit✔️✔️:             tmp = 9.23;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 13.14;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (tmp < 0.0) {
    // RDKit✔️✔️:         tmp = 28.5 - nNbrs[i] * 8.6 + nHs[i] * 1.5;
    // RDKit✔️✔️:         if (tmp < 0) {
    // RDKit✔️✔️:           tmp = 0.0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    let mut tmp = -1.0;
    match n_nbrs {
        1 => {
            if n_hs == 0 && chg == 0 && n_doub == 1 {
                tmp = 17.07;
            } else if n_hs == 1 && chg == 0 && n_sing == 1 {
                tmp = 20.23;
            } else if n_hs == 0 && chg == -1 && n_sing == 1 {
                tmp = 23.06;
            }
        }
        2 => {
            if n_hs == 0 && chg == 0 && n_sing == 2 && in3_ring {
                tmp = 12.53;
            } else if n_hs == 0 && chg == 0 && n_sing == 2 && !in3_ring {
                tmp = 9.23;
            } else if n_hs == 0 && chg == 0 && n_arom == 2 {
                tmp = 13.14;
            }
        }
        _ => {}
    }
    if tmp < 0.0 {
        tmp = 28.5 - f64::from(n_nbrs) * 8.6 + f64::from(n_hs) * 1.5;
        if tmp < 0.0 {
            tmp = 0.0;
        }
    }
    tmp
}

fn rdkit_tpsa_phosphorus_contrib(
    n_nbrs: i32,
    n_sing: i32,
    n_doub: i32,
    n_hs: i32,
    chg: i32,
) -> f64 {
    // RDKit✔️✔️:     } else if (includeSandP && atNum == 15) {
    // RDKit✔️✔️:       tmp = 0.0;
    // RDKit✔️✔️:       switch (nNbrs[i]) {
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 34.14;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3) {
    // RDKit✔️✔️:             tmp = 13.59;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 &&
    // RDKit✔️✔️:                      nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 23.47;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 9.81;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    let mut tmp = 0.0;
    match n_nbrs {
        2 => {
            if n_hs == 0 && chg == 0 && n_sing == 1 && n_doub == 1 {
                tmp = 34.14;
            }
        }
        3 => {
            if n_hs == 0 && chg == 0 && n_sing == 3 {
                tmp = 13.59;
            } else if n_hs == 1 && chg == 0 && n_sing == 2 && n_doub == 1 {
                tmp = 23.47;
            }
        }
        4 => {
            if n_hs == 0 && chg == 0 && n_sing == 3 && n_doub == 1 {
                tmp = 9.81;
            }
        }
        _ => {}
    }
    tmp
}

fn rdkit_tpsa_sulfur_contrib(
    n_nbrs: i32,
    n_sing: i32,
    n_doub: i32,
    n_arom: i32,
    n_hs: i32,
    chg: i32,
) -> f64 {
    // RDKit✔️✔️:     } else if (includeSandP && atNum == 16) {
    // RDKit✔️✔️:       tmp = 0.0;
    // RDKit✔️✔️:       switch (nNbrs[i]) {
    // RDKit✔️✔️:         case 1:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 32.09;
    // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 1) {
    // RDKit✔️✔️:             tmp = 38.80;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2) {
    // RDKit✔️✔️:             tmp = 25.30;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
    // RDKit✔️✔️:             tmp = 28.24;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nArom[i] == 2 && nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 21.70;
    // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 &&
    // RDKit✔️✔️:                      nDoub[i] == 1) {
    // RDKit✔️✔️:             tmp = 19.21;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && nDoub[i] == 2) {
    // RDKit✔️✔️:             tmp = 8.38;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    let mut tmp = 0.0;
    match n_nbrs {
        1 => {
            if n_hs == 0 && chg == 0 && n_doub == 1 {
                tmp = 32.09;
            } else if n_hs == 1 && chg == 0 && n_sing == 1 {
                tmp = 38.80;
            }
        }
        2 => {
            if n_hs == 0 && chg == 0 && n_sing == 2 {
                tmp = 25.30;
            } else if n_hs == 0 && chg == 0 && n_arom == 2 {
                tmp = 28.24;
            }
        }
        3 => {
            if n_hs == 0 && chg == 0 && n_arom == 2 && n_doub == 1 {
                tmp = 21.70;
            } else if n_hs == 0 && chg == 0 && n_sing == 2 && n_doub == 1 {
                tmp = 19.21;
            }
        }
        4 => {
            if n_hs == 0 && chg == 0 && n_sing == 2 && n_doub == 2 {
                tmp = 8.38;
            }
        }
        _ => {}
    }
    tmp
}

fn rdkit_tpsa_total_num_hs_without_neighbors(
    explicit_hydrogens: u8,
    implicit_hydrogens: i32,
) -> DescriptorResult<u32> {
    let total = i32::from(explicit_hydrogens) + implicit_hydrogens;
    u32::try_from(total).map_err(|_| DescriptorError::Unsupported {
        function: "calc_tpsa",
        rdkit_function: "Atom::getTotalNumHs",
    })
}

pub fn calc_num_aromatic_rings(mol: &Molecule) -> DescriptorResult<u32> {
    // RDKit✔️✔️: unsigned int calcNumAromaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let rings = symmetrize_sssr(mol).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aromatic_rings",
        rdkit_function: "Descriptors::calcNumAromaticRings/RingInfo::bondRings",
    })?;
    let mut res = 0_u32;
    for ring in rings.bond_rings() {
        res += 1;
        for bond_id in ring {
            if !mol.bonds()[bond_id.index()].is_aromatic() {
                res -= 1;
                break;
            }
        }
    }
    Ok(res)
}

pub fn calc_num_rotatable_bonds(
    mol: &Molecule,
    options: NumRotatableBondsOptions,
) -> DescriptorResult<u32> {
    // RDKit✔️✔️: #ifdef RDK_USE_STRICT_ROTOR_DEFINITION
    // RDKit✔️✔️: const NumRotatableBondsOptions DefaultStrictDefinition = Strict;
    // RDKit✔️✔️: #else
    // RDKit✔️✔️: const NumRotatableBondsOptions DefaultStrictDefinition = NonStrict;
    // RDKit✔️✔️: #endif
    //
    // RDKit✔️✔️: unsigned int calcNumRotatableBonds(const ROMol &mol,
    // RDKit✔️✔️:                                    NumRotatableBondsOptions strict) {
    // RDKit✔️✔️:   if (strict == Default) {
    // RDKit✔️✔️:     strict = DefaultStrictDefinition;
    // RDKit✔️✔️:   }
    //
    // RDKit✔️✔️:   if (strict == NonStrict) {
    // RDKit✔️✔️:     std::string pattern = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
    // RDKit✔️✔️:     pattern_flyweight m(pattern);
    // RDKit✔️✔️:     return m.get().countMatches(mol);
    // RDKit✔️✔️:   } else if (strict == Strict) {
    // RDKit✔️✔️:     std::string strict_pattern =
    // RDKit✔️✔️:         "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
    // RDKit✔️✔️:         "[CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
    // RDKit✔️✔️:         "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
    // RDKit✔️✔️:         "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
    // RDKit✔️✔️:         "CH3])[CH3])&!$([CH3])]";
    // RDKit✔️✔️:     pattern_flyweight m(strict_pattern);
    // RDKit✔️✔️:     return m.get().countMatches(mol);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     pattern_flyweight rotBonds_matcher("[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]");
    // RDKit✔️✔️:     pattern_flyweight nonRingAmides_matcher("[C&!R](=O)NC");
    // RDKit✔️✔️:     pattern_flyweight symRings_matcher(
    // RDKit✔️✔️:         "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])("
    // RDKit✔️✔️:         "a[!#1])a)]");
    // RDKit✔️✔️:     pattern_flyweight terminalTripleBonds_matcher("C#[#6,#7]");
    // RDKit✔️✔️:
    // RDKit✔️✔️:     std::vector<MatchVectType> matches;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int res = rotBonds_matcher.get().countMatches(mol);
    // RDKit✔️✔️:     if (!res) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     res -= symRings_matcher.get().countMatches(mol);
    // RDKit✔️✔️:     if (res < 0) {
    // RDKit✔️✔️:       res = 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     res -= terminalTripleBonds_matcher.get().countMatches(mol);
    // RDKit✔️✔️:     if (res < 0) {
    // RDKit✔️✔️:       res = 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
    // RDKit✔️✔️:     SubstructMatch(mol, *(nonRingAmides_matcher.get().getMatcher()), matches);
    // RDKit✔️✔️:     for (const auto &iv : matches) {
    // RDKit✔️✔️:       bool distinct = true;
    // RDKit✔️✔️:       for (const auto &mIt : iv) {
    // RDKit✔️✔️:         if (atomsSeen[mIt.second]) {
    // RDKit✔️✔️:           distinct = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atomsSeen.set(mIt.second);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (distinct && res > 0) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (res < 0) {
    // RDKit✔️✔️:       res = 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return static_cast<unsigned int>(res);
    // RDKit✔️✔️:   }
    if options == NumRotatableBondsOptions::StrictLinkages {
        return calc_num_rotatable_bonds_strict_linkages(mol);
    }
    let pattern = match options {
        NumRotatableBondsOptions::Default | NumRotatableBondsOptions::Strict => {
            RDKIT_ROTATABLE_BONDS_STRICT_SMARTS
        }
        NumRotatableBondsOptions::NonStrict => RDKIT_ROTATABLE_BONDS_NON_STRICT_SMARTS,
        NumRotatableBondsOptions::StrictLinkages => unreachable!("handled above"),
    };
    rdkit_count_smarts_matches("calc_num_rotatable_bonds", pattern)
        .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_rotatable_bonds"))
}

fn calc_num_rotatable_bonds_strict_linkages(mol: &Molecule) -> DescriptorResult<u32> {
    let mut res = rdkit_count_smarts_matches(
        "calc_num_rotatable_bonds",
        RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_BASE_SMARTS,
    )
    .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_rotatable_bonds"))?
        as i32;
    if res == 0 {
        return Ok(0);
    }

    res -= rdkit_count_smarts_matches(
        "calc_num_rotatable_bonds",
        RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_SYM_RINGS_SMARTS,
    )
    .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_rotatable_bonds"))?
        as i32;
    if res < 0 {
        res = 0;
    }

    res -= rdkit_count_smarts_matches(
        "calc_num_rotatable_bonds",
        RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_TERMINAL_TRIPLE_BONDS_SMARTS,
    )
    .and_then(|query| rdkit_count_query_matches(mol, &query, "calc_num_rotatable_bonds"))?
        as i32;
    if res < 0 {
        res = 0;
    }

    let amides = rdkit_count_smarts_matches(
        "calc_num_rotatable_bonds",
        RDKIT_ROTATABLE_BONDS_STRICT_LINKAGES_NON_RING_AMIDES_SMARTS,
    )?;
    let matches = crate::get_substruct_matches(mol, &amides);
    let mut atoms_seen = vec![false; mol.num_atoms()];
    for matched in matches {
        if rdkit_strict_linkages_amide_match_is_distinct(&mut atoms_seen, &matched) && res > 0 {
            res -= 1;
        }
    }

    if res < 0 {
        res = 0;
    }
    Ok(res as u32)
}

fn rdkit_strict_linkages_amide_match_is_distinct(
    atoms_seen: &mut [bool],
    matched: &SubstructMatchResult,
) -> bool {
    let mut distinct = true;
    for &atom_idx in &matched.atom_mapping {
        if atom_idx >= atoms_seen.len() {
            continue;
        }
        if atoms_seen[atom_idx] {
            distinct = false;
        }
        atoms_seen[atom_idx] = true;
    }
    distinct
}

#[derive(Debug, Clone, Copy)]
struct QedProperties {
    mw: f64,
    alogp: f64,
    hba: f64,
    hbd: f64,
    psa: f64,
    rotb: f64,
    arom: f64,
    alerts: f64,
}

impl QedProperties {
    fn values_in_rdkit_order(self) -> [f64; 8] {
        // RDKit✔️✔️: QEDproperties = namedtuple('QEDproperties', 'MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS')
        [
            self.mw,
            self.alogp,
            self.hba,
            self.hbd,
            self.psa,
            self.rotb,
            self.arom,
            self.alerts,
        ]
    }
}

#[derive(Debug, Clone, Copy)]
struct QedAdsParameter {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    dmax: f64,
}

const RDKIT_QED_WEIGHT_MEAN: QedProperties = QedProperties {
    // RDKit✔️✔️: QEDproperties = namedtuple('QEDproperties', 'MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS')
    // RDKit✔️✔️: WEIGHT_MEAN = QEDproperties(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)
    mw: 0.66,
    alogp: 0.46,
    hba: 0.05,
    hbd: 0.61,
    psa: 0.06,
    rotb: 0.65,
    arom: 0.48,
    alerts: 0.95,
};

const RDKIT_QED_ALIPHATIC_RINGS_SMARTS: &str = "[$([A;R][!a])]";
const RDKIT_QED_ACCEPTOR_SMARTS: &[&str] = &[
    // RDKit✔️✔️: AcceptorSmarts = [
    // RDKit✔️✔️:   '[oH0;X2]', '[OH1;X2;v2]', '[OH0;X2;v2]', '[OH0;X1;v2]', '[O-;X1]', '[SH0;X2;v2]', '[SH0;X1;v2]',
    // RDKit✔️✔️:   '[S-;X1]', '[nH0;X2]', '[NH0;X1;v3]', '[$([N;+0;X3;v3]);!$(N[C,S]=O)]'
    // RDKit✔️✔️: ]
    "[oH0;X2]",
    "[OH1;X2;v2]",
    "[OH0;X2;v2]",
    "[OH0;X1;v2]",
    "[O-;X1]",
    "[SH0;X2;v2]",
    "[SH0;X1;v2]",
    "[S-;X1]",
    "[nH0;X2]",
    "[NH0;X1;v3]",
    "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
];

const RDKIT_QED_STRUCTURAL_ALERT_SMARTS: &[&str] = &[
    // RDKit✔️✔️: StructuralAlertSmarts = [
    "*1[O,S,N]*1",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "[CX4][Cl,Br,I]",
    "[#6]S(=O)(=O)O[#6]",
    "[$([CH]),$(CC)]#CC(=O)[#6]",
    "[$([CH]),$(CC)]#CC(=O)O[#6]",
    "n[OH]",
    "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
    "C=C(C=O)C=O",
    "n1c([F,Cl,Br,I])cccc1",
    "[CH1](=O)",
    "[#8][#8]",
    "[C;!R]=[N;!R]",
    "[N!R]=[N!R]",
    "[#6](=O)[#6](=O)",
    "[#16][#16]",
    "[#7][NH2]",
    "C(=O)N[NH2]",
    "[#6]=S",
    "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    "C1(=[O,N])C=CC(=[O,N])C=C1",
    "C1(=[O,N])C(=[O,N])C=CC=C1",
    "a21aa3a(aa1aaaa2)aaaa3",
    "a31a(a2a(aa1)aaaa2)aaaa3",
    "a1aa2a3a(a1)A=AA=A3=AA=A2",
    "c1cc([NH2])ccc1",
    "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
    "I",
    "OS(=O)(=O)[O-]",
    "[N+](=O)[O-]",
    "C(=O)N[OH]",
    "C1NC(=O)NC(=O)1",
    "[SH]",
    "[S-]",
    "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
    "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
    "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
    "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "C#C",
    "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
    "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
    "[#6]=N[OH]",
    "[#6]=NOC=O",
    "[#6](=O)[CX4,CR0X3,O][#6](=O)",
    "c1ccc2c(c1)ccc(=O)o2",
    "[O+,o+,S+,s+]",
    "N=C=O",
    "[NX3,NX4][F,Cl,Br,I]",
    "c1ccccc1OC(=O)[#6]",
    "[CR0]=[CR0][CR0]=[CR0]",
    "[C+,c+,C-,c-]",
    "N=[N+]=[N-]",
    "C12C(NC(N1)=O)CSC2",
    "c1c([OH])c([OH,NH2,NH])ccc1",
    "P",
    "[N,O,S]C#N",
    "C=C=O",
    "[Si][F,Cl,Br,I]",
    "[SX2]O",
    "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
    "O1CCCCC1OC2CCC3CCCCC3C2",
    "N=[CR0][N,n,O,S]",
    "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
    "C=[C!r]C#N",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
    "[OH]c1ccc([OH,NH2,NH])cc1",
    "c1ccccc1OC(=O)O",
    "[SX2H0][N]",
    "c12ccccc1(SC(S)=N2)",
    "c12ccccc1(SC(=S)N2)",
    "c1nnnn1C=O",
    "s1c(S)nnc1NC=O",
    "S1C=CSC1=S",
    "C(=O)Onnn",
    "OS(=O)(=O)C(F)(F)F",
    "N#CC[OH]",
    "N#CC(=O)",
    "S(=O)(=O)C#N",
    "N[CH2]C#N",
    "C1(=O)NCC1",
    "S(=O)(=O)[O-,OH]",
    "NC[F,Cl,Br,I]",
    "C=[C!r]O",
    "[NX2+0]=[O+0]",
    "[OR0,NR0][OR0,NR0]",
    "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
    "[CX2R0][NX3R0]",
    "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
    "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
    "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
    "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
    "[*]=[N+]=[*]",
    "[SX3](=O)[O-,OH]",
    "N#N",
    "F.F.F.F",
    "[R0;D2][R0;D2][R0;D2][R0;D2]",
    "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
    "C=!@CC=[O,S]",
    "[#6,#8,#16][#6](=O)O[#6]",
    "c[C;R0](=[O,S])[#6]",
    "c[SX2][C;!R]",
    "C=C=C",
    "c1nc([F,Cl,Br,I,S])ncc1",
    "c1ncnc([F,Cl,Br,I,S])c1",
    "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
    "[#6]S(=O)(=O)c1ccc(cc1)F",
    "[15N]",
    "[13C]",
    "[18O]",
    "[34S]",
];

const RDKIT_QED_ADS_PARAMETERS: [QedAdsParameter; 8] = [
    // RDKit✔️✔️: adsParameters = {
    // RDKit✔️✔️:   'MW':
    QedAdsParameter {
        a: 2.817065973,
        b: 392.5754953,
        c: 290.7489764,
        d: 2.419764353,
        e: 49.22325677,
        f: 65.37051707,
        dmax: 104.9805561,
    },
    // RDKit✔️✔️:   'ALOGP':
    QedAdsParameter {
        a: 3.172690585,
        b: 137.8624751,
        c: 2.534937431,
        d: 4.581497897,
        e: 0.822739154,
        f: 0.576295591,
        dmax: 131.3186604,
    },
    // RDKit✔️✔️:   'HBA':
    QedAdsParameter {
        a: 2.948620388,
        b: 160.4605972,
        c: 3.615294657,
        d: 4.435986202,
        e: 0.290141953,
        f: 1.300669958,
        dmax: 148.7763046,
    },
    // RDKit✔️✔️:   'HBD':
    QedAdsParameter {
        a: 1.618662227,
        b: 1010.051101,
        c: 0.985094388,
        d: 0.000000001,
        e: 0.713820843,
        f: 0.920922555,
        dmax: 258.1632616,
    },
    // RDKit✔️✔️:   'PSA':
    QedAdsParameter {
        a: 1.876861559,
        b: 125.2232657,
        c: 62.90773554,
        d: 87.83366614,
        e: 12.01999824,
        f: 28.51324732,
        dmax: 104.5686167,
    },
    // RDKit✔️✔️:   'ROTB':
    QedAdsParameter {
        a: 0.010000000,
        b: 272.4121427,
        c: 2.558379970,
        d: 1.565547684,
        e: 1.271567166,
        f: 2.758063707,
        dmax: 105.4420403,
    },
    // RDKit✔️✔️:   'AROM':
    QedAdsParameter {
        a: 3.217788970,
        b: 957.7374108,
        c: 2.274627939,
        d: 0.000000001,
        e: 1.317690384,
        f: 0.375760881,
        dmax: 312.3372610,
    },
    // RDKit✔️✔️:   'ALERTS':
    QedAdsParameter {
        a: 0.010000000,
        b: 1199.094025,
        c: -0.09002883,
        d: 0.000000001,
        e: 0.185904477,
        f: 0.875193782,
        dmax: 417.7253140,
    },
];

fn rdkit_qed_ads(x: f64, p: QedAdsParameter) -> f64 {
    // RDKit✔️✔️: def ads(x, adsParameter):
    // RDKit✔️✔️:   """ ADS function """
    // RDKit✔️✔️:   p = adsParameter
    // RDKit✔️✔️:   exp1 = 1 + math.exp(-1 * (x - p.C + p.D / 2) / p.E)
    let exp1 = 1.0 + (-1.0 * (x - p.c + p.d / 2.0) / p.e).exp();
    // RDKit✔️✔️:   exp2 = 1 + math.exp(-1 * (x - p.C - p.D / 2) / p.F)
    let exp2 = 1.0 + (-1.0 * (x - p.c - p.d / 2.0) / p.f).exp();
    // RDKit✔️✔️:   dx = p.A + p.B / exp1 * (1 - 1 / exp2)
    let dx = p.a + p.b / exp1 * (1.0 - 1.0 / exp2);
    // RDKit✔️✔️:   return dx / p.DMAX
    dx / p.dmax
}

fn rdkit_qed_properties(mol: &Molecule) -> DescriptorResult<QedProperties> {
    // RDKit✔️✔️: def properties(mol):
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   Calculates the properties that are required to calculate the QED descriptor.
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   if mol is None:
    // RDKit✔️✔️:     raise ValueError('You need to provide a mol argument.')
    // COSMolKit modeled input state: Rust core APIs receive `&Molecule`, so
    // Python's `None` input is unrepresentable at this boundary. The compile
    // time API-boundary test below fixes this as an excluded input state.
    // RDKit✔️✔️:   mol = Chem.RemoveHs(mol)
    let mol_no_h = mol
        .without_hydrogens()
        .map_err(|_| DescriptorError::Unsupported {
            function: "calc_qed",
            rdkit_function: "Chem.RemoveHs",
        })?;
    // RDKit✔️✔️:   qedProperties = QEDproperties(
    // RDKit✔️✔️:     MW=rdmd._CalcMolWt(mol),
    let mw = calc_mol_wt(&mol_no_h, false)?;
    // RDKit✔️✔️:     ALOGP=Crippen.MolLogP(mol),
    let alogp = calc_crippen_descriptors(&mol_no_h, true, false)?.logp;
    // RDKit✔️✔️:     HBA=sum(
    // RDKit✔️✔️:       len(mol.GetSubstructMatches(pattern)) for pattern in Acceptors
    // RDKit✔️✔️:       if mol.HasSubstructMatch(pattern)),
    let hba = rdkit_qed_hba(&mol_no_h)?;
    // RDKit✔️✔️:     HBD=rdmd.CalcNumHBD(mol),
    let hbd = calc_num_hbd(&mol_no_h)?;
    // RDKit✔️✔️:     PSA=MolSurf.TPSA(mol),
    let psa = calc_tpsa(&mol_no_h, false, false)?;
    // RDKit✔️✔️:     ROTB=rdmd.CalcNumRotatableBonds(mol, rdmd.NumRotatableBondsOptions.Strict),
    let rotb = calc_num_rotatable_bonds(&mol_no_h, NumRotatableBondsOptions::Strict)?;
    // RDKit✔️✔️:     AROM=len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), AliphaticRings))),
    let arom = rdkit_qed_arom(&mol_no_h)?;
    // RDKit✔️✔️:     ALERTS=sum(1 for alert in StructuralAlerts if mol.HasSubstructMatch(alert)),
    let alerts = rdkit_qed_alerts(&mol_no_h)?;
    // RDKit✔️✔️:   )
    // RDKit✔️✔️:   return qedProperties
    Ok(QedProperties {
        mw,
        alogp,
        hba: f64::from(hba),
        hbd: f64::from(hbd),
        psa,
        rotb: f64::from(rotb),
        arom: f64::from(arom),
        alerts: f64::from(alerts),
    })
}

fn rdkit_qed_hba(mol: &Molecule) -> DescriptorResult<u32> {
    let mut count = 0_u32;
    for &pattern in RDKIT_QED_ACCEPTOR_SMARTS {
        count = count
            .checked_add(rdkit_qed_substruct_match_count(mol, pattern)?)
            .ok_or(DescriptorError::Unsupported {
                function: "calc_qed",
                rdkit_function: "QED.properties HBA sum",
            })?;
    }
    Ok(count)
}

fn rdkit_qed_alerts(mol: &Molecule) -> DescriptorResult<u32> {
    let mut count = 0_u32;
    for &pattern in RDKIT_QED_STRUCTURAL_ALERT_SMARTS {
        let query = rdkit_count_smarts_matches("calc_qed", pattern)?;
        if crate::has_substruct_match(mol, &query) {
            count += 1;
        }
    }
    Ok(count)
}

fn rdkit_qed_substruct_match_count(mol: &Molecule, pattern: &str) -> DescriptorResult<u32> {
    let query = rdkit_count_smarts_matches("calc_qed", pattern)?;
    if !crate::has_substruct_match(mol, &query) {
        return Ok(0);
    }
    let matches = crate::get_substruct_matches(mol, &query);
    u32::try_from(matches.len()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_qed",
        rdkit_function: "ROMol.GetSubstructMatches",
    })
}

fn rdkit_qed_arom(mol: &Molecule) -> DescriptorResult<u32> {
    let aliphatic_rings = rdkit_count_smarts_matches("calc_qed", RDKIT_QED_ALIPHATIC_RINGS_SMARTS)?;
    let mol_without_aliphatic_rings =
        rdkit_qed_delete_substructs(mol, &aliphatic_rings, false, false)?;
    let rings = symmetrize_sssr(&mol_without_aliphatic_rings).map_err(|_| {
        DescriptorError::Unsupported {
            function: "calc_qed",
            rdkit_function: "Chem.GetSSSR",
        }
    })?;
    u32::try_from(rings.atom_rings().len()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_qed",
        rdkit_function: "Chem.GetSSSR",
    })
}

fn rdkit_qed_delete_substructs(
    mol: &Molecule,
    query: &Molecule,
    only_frags: bool,
    use_chirality: bool,
) -> DescriptorResult<Molecule> {
    // RDKit✔️✔️: ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query, bool onlyFrags,
    // RDKit✔️✔️:                         bool useChirality) {
    // RDKit✔️✔️:   auto *res = new RWMol(mol, false);
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   // do the substructure matching and get the atoms that match the query
    // RDKit✔️✔️:   const bool uniquify = true;
    // RDKit✔️✔️:   const bool recursionPossible = true;
    // RDKit✔️✔️:   SubstructMatch(*res, query, fgpMatches, uniquify, recursionPossible,
    // RDKit✔️✔️:                  useChirality);
    let mut params = SubstructMatchParams::default();
    params.use_chirality = use_chirality;
    let fgp_matches =
        crate::try_get_substruct_matches_with_params(mol, query, &params).map_err(|_| {
            DescriptorError::Unsupported {
                function: "Chem.DeleteSubstructs",
                rdkit_function: "SubstructMatch(..., useChirality=true)",
            }
        })?;

    // RDKit✔️✔️:   // if didn't find any matches nothing to be done here
    // RDKit✔️✔️:   // simply return a copy of the molecule
    // RDKit✔️✔️:   if (fgpMatches.empty()) {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    if fgp_matches.is_empty() {
        return Ok(mol.clone());
    }

    // RDKit✔️✔️:   // all matches on the molecule - list of list of atom ids
    // RDKit✔️✔️:   VECT_INT_VECT matches;
    // RDKit✔️✔️:   matches.reserve(fgpMatches.size());
    // RDKit✔️✔️:   for (const auto &mati : fgpMatches) {
    // RDKit✔️✔️:     INT_VECT match;  // each match onto the molecule - list of atoms ids
    // RDKit✔️✔️:     match.reserve(mati.size());
    // RDKit✔️✔️:     for (const auto &mi : mati) {
    // RDKit✔️✔️:       match.push_back(mi.second);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     matches.push_back(std::move(match));
    // RDKit✔️✔️:   }
    //
    // RDKit✔️✔️:   // now loop over the list of matches and check if we can delete any of them
    // RDKit✔️✔️:   INT_VECT delList;
    // RDKit✔️✔️:   if (onlyFrags) {
    // RDKit✔️✔️:     VECT_INT_VECT frags;
    // RDKit✔️✔️:     MolOps::getMolFrags(*res, frags);
    // RDKit✔️✔️:     for (auto &fi : frags) {
    // RDKit✔️✔️:       std::sort(fi.begin(), fi.end());
    // RDKit✔️✔️:       for (auto &mxi : matches) {
    // RDKit✔️✔️:         std::sort(mxi.begin(), mxi.end());
    // RDKit✔️✔️:         if (fi == mxi) {
    // RDKit✔️✔️:           INT_VECT tmp;
    // RDKit✔️✔️:           Union(mxi, delList, tmp);
    // RDKit✔️✔️:           delList = tmp;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }  // end of if we found a matching fragment
    // RDKit✔️✔️:       }  // end of loop over matches
    // RDKit✔️✔️:     }  // end of loop over fragments
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // in this case we want to delete any matches we find
    // RDKit✔️✔️:     // simply loop over the matches and collect the atoms that need to
    // RDKit✔️✔️:     // be removed
    // RDKit✔️✔️:     for (const auto &mxi : matches) {
    // RDKit✔️✔️:       INT_VECT tmp;
    // RDKit✔️✔️:       Union(mxi, delList, tmp);
    // RDKit✔️✔️:       delList = tmp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut del_atoms = BTreeSet::<usize>::new();
    if only_frags {
        let fragments = rdkit_delete_substructs_fragments(mol);
        let mut matches = fgp_matches
            .iter()
            .map(|matched| matched.atom_mapping.iter().copied().collect::<Vec<_>>())
            .collect::<Vec<_>>();
        for mut fragment in fragments {
            fragment.sort_unstable();
            for matched in &mut matches {
                matched.sort_unstable();
                if fragment == *matched {
                    del_atoms.extend(matched.iter().copied());
                    break;
                }
            }
        }
    } else {
        for matched in &fgp_matches {
            for &atom_idx in &matched.atom_mapping {
                del_atoms.insert(atom_idx);
            }
        }
    }
    let del_list = del_atoms.into_iter().map(AtomId::new).collect::<Vec<_>>();

    // RDKit✔️✔️:   if (delList.empty()) {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    if del_list.is_empty() {
        return Ok(mol.clone());
    }

    // RDKit✔️✔️:   // now loop over the union list and delete the atoms
    // RDKit✔️✔️:   res->beginBatchEdit();
    // RDKit✔️✔️:   boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   for (auto idx : delList) {
    // RDKit✔️✔️:     removedAtoms.set(idx);
    // RDKit✔️✔️:     res->removeAtom(idx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res->commitBatchEdit();
    let mut builder = mol.to_builder();
    builder.remove_atoms_for_construction(&del_list);

    // RDKit✔️✔️:   details::updateSubMolConfs(mol, *res, removedAtoms);
    // RDKit✔️✔️:   res->clearComputedProps(true);
    // RDKit✔️✔️:   // update our properties, but allow unhappiness:
    // RDKit✔️✔️:   res->updatePropertyCache(false);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    builder.build().map_err(|_| DescriptorError::Unsupported {
        function: "calc_qed",
        rdkit_function: "Chem.DeleteSubstructs",
    })
}

fn rdkit_delete_substructs_fragments(mol: &Molecule) -> Vec<Vec<usize>> {
    // RDKit✔️✔️: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit✔️✔️:   unsigned int natms = mol.getNumAtoms();
    // RDKit✔️✔️:   mapping.resize(natms);
    // RDKit✔️✔️:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit✔️✔️:                : 0;
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
    // RDKit✔️✔️:   frags.clear();
    // RDKit✔️✔️:   INT_VECT mapping;
    // RDKit✔️✔️:   getMolFrags(mol, mapping);
    let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    let mut mapping = vec![usize::MAX; mol.num_atoms()];
    let mut component_count = 0usize;
    for start in 0..mol.num_atoms() {
        if mapping[start] != usize::MAX {
            continue;
        }
        let mut queue = VecDeque::new();
        queue.push_back(start);
        mapping[start] = component_count;
        while let Some(atom) = queue.pop_front() {
            for neighbor in adjacency.neighbors_of(atom) {
                if mapping[neighbor.atom_index] == usize::MAX {
                    mapping[neighbor.atom_index] = component_count;
                    queue.push_back(neighbor.atom_index);
                }
            }
        }
        component_count += 1;
    }
    let mut components = BTreeMap::<usize, Vec<usize>>::new();
    // RDKit✔️✔️:   INT_INT_VECT_MAP comMap;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️✔️:     int mi = mapping[i];
    // RDKit✔️✔️:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️✔️:       INT_VECT comp;
    // RDKit✔️✔️:       comMap[mi] = comp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     comMap[mi].push_back(i);
    for (atom, component) in mapping.into_iter().enumerate() {
        components.entry(component).or_default().push(atom);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️✔️:     frags.push_back((*mci).second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return rdcast<unsigned int>(frags.size());
    // RDKit✔️✔️: }
    components.into_values().collect()
}

pub fn calc_qed(mol: &Molecule) -> DescriptorResult<f64> {
    // RDKit✔️✔️: def qed(mol, w=WEIGHT_MEAN, qedProperties=None):
    // RDKit✔️✔️:   """ Calculate the weighted sum of ADS mapped properties
    // RDKit✔️✔️:
    // RDKit✔️✔️:   some examples from the QED paper, reference values from Peter G's original implementation
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.253...
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('CNC(=NCCSCc1nc[nH]c1C)NC#N')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.234...
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.234...
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   if qedProperties is None:
    // RDKit✔️✔️:     qedProperties = properties(mol)
    let qed_properties = rdkit_qed_properties(mol)?;
    // RDKit✔️✔️:   d = [ads(pi, adsParameters[name]) for name, pi in qedProperties._asdict().items()]
    let values = qed_properties.values_in_rdkit_order();
    let d = values
        .into_iter()
        .zip(RDKIT_QED_ADS_PARAMETERS)
        .map(|(pi, parameter)| rdkit_qed_ads(pi, parameter))
        .collect::<Vec<_>>();
    // RDKit✔️✔️:   t = sum(wi * math.log(di) for wi, di in zip(w, d))
    let weights = RDKIT_QED_WEIGHT_MEAN.values_in_rdkit_order();
    let t = weights
        .iter()
        .copied()
        .zip(d.iter().copied())
        .map(|(wi, di)| wi * di.ln())
        .sum::<f64>();
    // RDKit✔️✔️:   return math.exp(t / sum(w))
    Ok((t / weights.iter().copied().sum::<f64>()).exp())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct FormulaKey {
    isotope: u32,
    symbol: &'static str,
}

fn hill_compare_ordering(v1: &FormulaKey, v2: &FormulaKey) -> Ordering {
    // RDKit✔️✔️: bool HillCompare(const std::pair<unsigned int, std::string> &v1,
    // RDKit✔️✔️:                  const std::pair<unsigned int, std::string> &v2) {
    // RDKit✔️✔️:   bool nCompare = (v1.first < v2.first);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // special cases: Cs, Hs, Ds, and Ts go at the beginning
    // RDKit✔️✔️:   if (v1.second == "C") {
    // RDKit✔️✔️:     if (v2.second != "C") {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return nCompare;
    // RDKit✔️✔️:   } else if (v2.second == "C") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (v1.second == "H") {
    // RDKit✔️✔️:     if (v2.second != "H") {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return nCompare;
    // RDKit✔️✔️:   } else if (v2.second == "H") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (v1.second == "D") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (v2.second == "D") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (v1.second == "T") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (v2.second == "T") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // finally, just compare the symbols and the isotopes:
    // RDKit✔️✔️:   if (v1 != v2) {
    // RDKit✔️✔️:     return v1 < v2;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return nCompare;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if v1.symbol == "C" {
        if v2.symbol != "C" {
            return Ordering::Less;
        }
        return v1.isotope.cmp(&v2.isotope);
    } else if v2.symbol == "C" {
        return Ordering::Greater;
    }
    if v1.symbol == "H" {
        if v2.symbol != "H" {
            return Ordering::Less;
        }
        return v1.isotope.cmp(&v2.isotope);
    } else if v2.symbol == "H" {
        return Ordering::Greater;
    }
    if v1.symbol == "D" {
        return Ordering::Less;
    } else if v2.symbol == "D" {
        return Ordering::Greater;
    }
    if v1.symbol == "T" {
        return Ordering::Less;
    } else if v2.symbol == "T" {
        return Ordering::Greater;
    }
    v1.cmp(v2)
}

fn rdkit_total_num_hs_without_neighbors(
    explicit_hydrogens: u8,
    implicit_hydrogens: i32,
) -> DescriptorResult<u32> {
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let total = i32::from(explicit_hydrogens) + implicit_hydrogens;
    u32::try_from(total).map_err(|_| DescriptorError::Unsupported {
        function: "calc_mol_formula",
        rdkit_function: "Atom::getTotalNumHs",
    })
}

fn rdkit_count_smarts_matches(function: &'static str, pattern: &str) -> DescriptorResult<Molecule> {
    // RDKit✔️✔️: ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    // RDKit✔️✔️:   m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    // RDKit✔️✔️:   RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit✔️✔️:   m_matcher = p;
    // RDKit✔️✔️:   POSTCONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️: };
    crate::chemistry::forcefield::crystalff::build_crystalff_query_molecule(pattern).map_err(|_| {
        DescriptorError::Unsupported {
            function,
            rdkit_function: "RDKit::SmartsToMol",
        }
    })
}

fn rdkit_count_query_matches(
    mol: &Molecule,
    query: &Molecule,
    function: &'static str,
) -> DescriptorResult<u32> {
    // RDKit✔️✔️: unsigned int countMatches(const RDKit::ROMol &mol) const {
    // RDKit✔️✔️:   PRECONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:   std::vector<RDKit::MatchVectType> matches;
    // RDKit✔️✔️:   if (m_needCopies) {
    // RDKit✔️✔️:     const RDKit::ROMol nm(*(m_matcher), true);
    // RDKit✔️✔️:     RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     const RDKit::ROMol &nm = *m_matcher;
    // RDKit✔️✔️:     RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return matches.size();
    // RDKit✔️✔️: }
    //
    // RDKit✔️✔️: bool recursionPossible = true;
    // RDKit✔️✔️: bool uniquify = true;
    // RDKit✔️✔️: unsigned int maxMatches = 1000;
    let matches = crate::get_substruct_matches(mol, query);
    u32::try_from(matches.len()).map_err(|_| DescriptorError::Unsupported {
        function,
        rdkit_function: "RDKit::SubstructMatch",
    })
}

fn unsupported<T>(function: &'static str, rdkit_function: &'static str) -> DescriptorResult<T> {
    Err(DescriptorError::Unsupported {
        function,
        rdkit_function,
    })
}

#[cfg(test)]
mod qed_tests {
    use super::*;
    use crate::SmilesWriteParams;
    use serde::Deserialize;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    #[derive(Debug, Deserialize)]
    struct DeleteSubstructsRecord {
        case: String,
        smiles: String,
        smarts: String,
        only_frags: bool,
        use_chirality: bool,
        rdkit_ok: bool,
        result_smiles: Option<String>,
        num_atoms: Option<usize>,
        num_bonds: Option<usize>,
        error: Option<String>,
    }

    fn load_delete_substructs_golden() -> Vec<DeleteSubstructsRecord> {
        let path = crate::test_data::golden_path("delete_substructs_onlyfrags_chirality.jsonl");
        assert!(
            path.exists(),
            "missing RDKit DeleteSubstructs golden: {}",
            path.display()
        );
        let file = File::open(&path).unwrap_or_else(|err| {
            panic!("failed to open {}: {err}", path.display());
        });
        BufReader::new(file)
            .lines()
            .enumerate()
            .map(|(idx, line)| {
                let line = line.unwrap_or_else(|err| {
                    panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
                });
                serde_json::from_str(&line).unwrap_or_else(|err| {
                    panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
                })
            })
            .collect()
    }

    fn structural_alert_hits(smiles: &str) -> Vec<usize> {
        let mol = Molecule::from_smiles(smiles).expect("test molecule must parse");
        let mol = mol
            .without_hydrogens()
            .expect("test molecule hydrogen removal must succeed");
        RDKIT_QED_STRUCTURAL_ALERT_SMARTS
            .iter()
            .enumerate()
            .filter_map(|(idx, pattern)| {
                let query =
                    rdkit_count_smarts_matches("calc_qed", pattern).expect("alert SMARTS parses");
                crate::has_substruct_match(&mol, &query).then_some(idx)
            })
            .collect()
    }

    #[test]
    fn qed_structural_alert_hits_match_rdkit_simple_smoke_rows() {
        assert_eq!(structural_alert_hits("C=C"), vec![19]);
        assert_eq!(structural_alert_hits("N#C"), Vec::<usize>::new());
        assert_eq!(structural_alert_hits("c1ccsc1"), Vec::<usize>::new());
        assert_eq!(structural_alert_hits("c1cc[se]c1"), vec![26]);
        assert_eq!(structural_alert_hits("CC(=O)O"), Vec::<usize>::new());
        assert_eq!(structural_alert_hits("F[C@@H]1O[C@H](Cl)S1"), vec![2]);
        assert_eq!(
            structural_alert_hits(
                "NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](CC(=O)N)NC(=O)[C@@H](N2)CCC(=O)N)CSCCCC(=O)N(C)[C@H](C(=O)N[C@H](C2=O)[C@@H](C)CC)Cc3ccc(O)cc3"
            ),
            Vec::<usize>::new()
        );
        assert_eq!(
            structural_alert_hits("O=[S](=O)(CCC(=O)N/C1=C/C=C(/NC(=O)C)C=C1)C2=CC=CC3=NON=C23"),
            Vec::<usize>::new()
        );
    }

    #[test]
    fn qed_properties_none_input_is_unrepresentable_in_rust_core_api() {
        let _: fn(&Molecule) -> DescriptorResult<f64> = calc_qed;
        let _: fn(&Molecule) -> DescriptorResult<QedProperties> = rdkit_qed_properties;
    }

    #[test]
    fn delete_substructs_golden_covers_required_parameter_branches() {
        let records = load_delete_substructs_golden();
        assert!(
            records.iter().any(|record| !record.only_frags),
            "DeleteSubstructs golden must cover onlyFrags=false"
        );
        assert!(
            records.iter().any(|record| record.only_frags),
            "DeleteSubstructs golden must cover onlyFrags=true"
        );
        assert!(
            records.iter().any(|record| record.use_chirality),
            "DeleteSubstructs golden must cover useChirality=true"
        );
        assert!(
            records.iter().any(|record| !record.use_chirality),
            "DeleteSubstructs golden must cover useChirality=false"
        );
    }

    #[test]
    fn delete_substructs_matches_rdkit_golden_for_only_frags_matrix() {
        let records = load_delete_substructs_golden();
        let mut failures = Vec::<String>::new();
        for (row_idx, record) in records.iter().enumerate() {
            let context = format!(
                "row {} case {} smiles={} smarts={} onlyFrags={} useChirality={}",
                row_idx + 1,
                record.case,
                record.smiles,
                record.smarts,
                record.only_frags,
                record.use_chirality
            );
            if !record.rdkit_ok {
                assert!(
                    record.error.is_some(),
                    "RDKit-not-ok DeleteSubstructs record has no error in {context}"
                );
                continue;
            }
            let expected_smiles = record.result_smiles.as_ref().unwrap_or_else(|| {
                panic!("RDKit-ok DeleteSubstructs record missing result_smiles in {context}")
            });
            let expected_atoms = record.num_atoms.unwrap_or_else(|| {
                panic!("RDKit-ok DeleteSubstructs record missing num_atoms in {context}")
            });
            let expected_bonds = record.num_bonds.unwrap_or_else(|| {
                panic!("RDKit-ok DeleteSubstructs record missing num_bonds in {context}")
            });
            let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
                panic!("COSMolKit failed to parse molecule in {context}: {err}")
            });
            let query = rdkit_count_smarts_matches("Chem.DeleteSubstructs", &record.smarts)
                .unwrap_or_else(|err| {
                    panic!("COSMolKit failed to parse SMARTS in {context}: {err}")
                });
            let actual =
                rdkit_qed_delete_substructs(&mol, &query, record.only_frags, record.use_chirality);
            match actual {
                Ok(actual) => {
                    let actual_smiles = actual
                        .to_smiles_with_params(&SmilesWriteParams {
                            do_isomeric_smiles: true,
                            canonical: true,
                            ..Default::default()
                        })
                        .unwrap_or_else(|err| {
                            panic!(
                                "COSMolKit failed to write DeleteSubstructs result in {context}: {err}"
                            )
                        });
                    if actual_smiles != *expected_smiles {
                        failures.push(format!(
                            "{context}: result_smiles mismatch actual={actual_smiles:?} expected={expected_smiles:?}"
                        ));
                    }
                    if actual.num_atoms() != expected_atoms {
                        failures.push(format!(
                            "{context}: num_atoms mismatch actual={} expected={expected_atoms}",
                            actual.num_atoms()
                        ));
                    }
                    if actual.num_bonds() != expected_bonds {
                        failures.push(format!(
                            "{context}: num_bonds mismatch actual={} expected={expected_bonds}",
                            actual.num_bonds()
                        ));
                    }
                }
                Err(err) => failures.push(format!("{context}: failed closed: {err}")),
            }
        }
        if !failures.is_empty() {
            let sample = failures
                .iter()
                .take(24)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n");
            panic!(
                "DeleteSubstructs parity failed with {} failures\nfirst failures:\n{sample}",
                failures.len()
            );
        }
    }
}
