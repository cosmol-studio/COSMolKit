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
}
