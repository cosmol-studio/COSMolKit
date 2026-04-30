/// Bond order for chemical edges.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Dative,
    Null,
}

/// RDKit-style directional single-bond marker used for SMILES cis/trans stereo.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondDirection {
    None,
    EndUpRight,
    EndDownRight,
}

/// RDKit-style double-bond stereo assignment.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondStereo {
    None,
    Any,
    Cis,
    Trans,
}

/// Bond record in a molecule graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    /// 0-based index in molecule bond table.
    pub index: usize,
    /// 0-based begin atom index.
    pub begin_atom: usize,
    /// 0-based end atom index.
    pub end_atom: usize,
    /// Bond order annotation.
    pub order: BondOrder,
    /// RDKit-style aromatic flag, independent from bond order.
    pub is_aromatic: bool,
    /// Directional single-bond marker used to assign double-bond stereo.
    pub direction: BondDirection,
    /// Double-bond stereo assignment after RDKit-like stereo perception.
    pub stereo: BondStereo,
    /// Controlling atom pair for double-bond stereo, in RDKit stereo atom order.
    pub stereo_atoms: Vec<usize>,
}
