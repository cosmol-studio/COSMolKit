#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ChiralTag {
    Unspecified,
    TetrahedralCw,
    TetrahedralCcw,
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
    /// Parsed tetrahedral chirality tag from SMILES (`@`/`@@`) when present.
    pub chiral_tag: ChiralTag,
    /// Optional isotope label from bracket SMILES.
    pub isotope: Option<u16>,
    /// Optional atom-map number from bracket SMILES.
    pub atom_map_num: Option<u32>,
}
