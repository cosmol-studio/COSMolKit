use std::collections::BTreeMap;

use crate::query::{BondQueryPredicate, QueryNode};

/// Bond order for chemical edges.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Dative,
    Hydrogen,
    Null,
}

macro_rules! impl_enum_metadata {
    (
        $type:ty,
        code $code_fn:ident,
        name $name_fn:ident,
        parse $parse_fn:ident,
        members [$($variant:path => ($code:literal, $name:literal, [$($alias:literal),* $(,)?])),+ $(,)?]
    ) => {
        impl $type {
            #[must_use]
            pub const fn $code_fn(self) -> i64 {
                match self {
                    $($variant => $code,)+
                }
            }

            #[must_use]
            pub const fn $name_fn(self) -> &'static str {
                match self {
                    $($variant => $name,)+
                }
            }

            #[must_use]
            pub fn $parse_fn(name: &str) -> Option<Self> {
                match name {
                    $($name $(| $alias)* => Some($variant),)+
                    _ => None,
                }
            }
        }
    };
}

impl_enum_metadata!(
    BondOrder,
    code python_code,
    name rdkit_name,
    parse from_rdkit_name,
    members [
        BondOrder::Null => (0, "UNSPECIFIED", ["ZERO"]),
        BondOrder::Single => (1, "SINGLE", []),
        BondOrder::Double => (2, "DOUBLE", []),
        BondOrder::Triple => (3, "TRIPLE", []),
        BondOrder::Quadruple => (4, "QUADRUPLE", []),
        BondOrder::Aromatic => (5, "AROMATIC", []),
        BondOrder::Dative => (6, "DATIVE", ["DATIVEL", "DATIVER"]),
        BondOrder::Hydrogen => (7, "HYDROGEN", []),
    ]
);

/// RDKit-style directional single-bond marker used for SMILES cis/trans stereo.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondDirection {
    None,
    EndUpRight,
    EndDownRight,
    Unknown,
}

impl_enum_metadata!(
    BondDirection,
    code python_code,
    name rdkit_name,
    parse from_rdkit_name,
    members [
        BondDirection::None => (0, "NONE", []),
        BondDirection::EndUpRight => (1, "ENDUPRIGHT", []),
        BondDirection::EndDownRight => (2, "ENDDOWNRIGHT", []),
        BondDirection::Unknown => (3, "UNKNOWN", []),
    ]
);

/// RDKit-style double-bond stereo assignment.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BondStereo {
    None,
    Any,
    Cis,
    Trans,
}

impl_enum_metadata!(
    BondStereo,
    code python_code,
    name rdkit_name,
    parse from_rdkit_name,
    members [
        BondStereo::None => (0, "STEREONONE", []),
        BondStereo::Any => (1, "STEREOANY", []),
        BondStereo::Cis => (2, "STEREOCIS", ["STEREOZ"]),
        BondStereo::Trans => (3, "STEREOTRANS", ["STEREOE"]),
    ]
);

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
    /// MDL molfile query-bond code for query bonds representable by a single
    /// V2000/V3000 bond type field.
    pub molfile_query_bond_code: Option<u8>,
    /// Preserved Molfile/SDF bond properties.
    pub props: BTreeMap<String, String>,
    /// RDKit-like query bond AST when the bond originates from a query Molfile.
    pub query: Option<QueryNode<BondQueryPredicate>>,
}

impl Bond {
    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }
}
