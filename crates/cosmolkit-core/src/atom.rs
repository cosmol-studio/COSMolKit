use std::collections::BTreeMap;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ChiralTag {
    Unspecified,
    TetrahedralCw,
    TetrahedralCcw,
    TrigonalBipyramidal,
}

macro_rules! impl_chiral_tag_metadata {
    ($($variant:path => ($code:literal, $name:literal, $from_rdkit:literal)),+ $(,)?) => {
        impl ChiralTag {
            #[must_use]
            pub const fn python_code(self) -> i64 {
                match self {
                    $($variant => $code,)+
                }
            }

            #[must_use]
            pub const fn rdkit_name(self) -> &'static str {
                match self {
                    $($variant => $name,)+
                }
            }

            #[must_use]
            pub fn from_rdkit_name(name: &str) -> Option<Self> {
                match name {
                    $($from_rdkit => Some($variant),)+
                    _ => None,
                }
            }
        }
    };
}

impl_chiral_tag_metadata! {
    ChiralTag::Unspecified => (0, "CHI_UNSPECIFIED", "CHI_UNSPECIFIED"),
    ChiralTag::TetrahedralCw => (1, "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CW"),
    ChiralTag::TetrahedralCcw => (2, "CHI_TETRAHEDRAL_CCW", "CHI_TETRAHEDRAL_CCW"),
    ChiralTag::TrigonalBipyramidal => (
        3,
        "CHI_TRIGONALBIPYRAMIDAL",
        "CHI_TRIGONALBIPYRAMIDAL"
    ),
}

/// Atom record in a molecule graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    /// 0-based index in molecule atom table.
    pub index: usize,
    /// Atomic number, e.g. 6 for carbon.
    pub atomic_num: u8,
    /// Aromatic flag as assigned during parsing.
    pub is_aromatic: bool,
    /// Formal charge in electron units.
    pub formal_charge: i8,
    /// Explicit hydrogen count provided in bracket SMILES.
    pub explicit_hydrogens: u8,
    /// RDKit-style noImplicit flag (true for bracket atoms in SMILES parser).
    pub no_implicit: bool,
    /// Radical electron count cache.
    pub num_radical_electrons: u8,
    /// Parsed or perceived RDKit chirality tag when present.
    pub chiral_tag: ChiralTag,
    /// Optional isotope label from bracket SMILES.
    pub isotope: Option<u16>,
    /// Optional atom-map number from bracket SMILES.
    pub atom_map_num: Option<u32>,
    /// Preserved Molfile/SDF atom properties.
    pub props: BTreeMap<String, String>,
    /// Cached RDKit `_CIPRank` atom property assigned by legacy stereochemistry.
    #[doc(hidden)]
    pub rdkit_cip_rank: Option<i64>,
}

impl Atom {
    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub fn prop_f64(&self, key: &str) -> Option<f64> {
        self.prop(key)?.parse().ok()
    }
}
