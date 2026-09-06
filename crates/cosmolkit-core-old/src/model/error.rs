use crate::{AtomId, BondId, BondStereo, SubstanceGroupId};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MoleculeBuildError {
    #[error("bond endpoint out of range: {begin}-{end}, atom_count={atom_count}")]
    BondEndpointOutOfRange {
        begin: AtomId,
        end: AtomId,
        atom_count: usize,
    },
    #[error("bond stereo {stereo:?} requires two stereo atom references")]
    BondStereoAtomsRequired { stereo: BondStereo },
    #[error("self-loop bond is not allowed: atom {atom}")]
    SelfLoopBond { atom: AtomId },
    #[error("2D coordinate row count mismatch: rows={rows}, atom_count={atom_count}")]
    CoordinateRowCount { rows: usize, atom_count: usize },
    #[error("3D conformer row count mismatch: rows={rows}, atom_count={atom_count}")]
    ConformerRowCount { rows: usize, atom_count: usize },
    #[error("3D conformer index out of range: index={index}, conformer_count={conformer_count}")]
    ConformerIndexOutOfRange {
        index: usize,
        conformer_count: usize,
    },
    #[error("substance group atom out of range: {atom}, atom_count={atom_count}")]
    SubstanceGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("substance group bond out of range: {bond}, bond_count={bond_count}")]
    SubstanceGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("substance group parent out of range: {parent:?}")]
    SubstanceGroupParentOutOfRange { parent: SubstanceGroupId },
    #[error("invalid molecule state produced during build: {message}")]
    InvalidMoleculeState { message: String },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum InvariantError {
    #[error("model value validation failed: {message}")]
    ModelValidation { message: String },
    #[error("atom table index mismatch: row {row} contains atom id {actual}")]
    AtomIndexMismatch { row: usize, actual: AtomId },
    #[error("bond table index mismatch: row {row} contains bond id {actual}")]
    BondIndexMismatch { row: usize, actual: BondId },
    #[error(
        "substance group table index mismatch: row {row} contains substance group id {actual:?}"
    )]
    SubstanceGroupIndexMismatch {
        row: usize,
        actual: SubstanceGroupId,
    },
    #[error("invalid bond endpoint on bond {bond}: {begin}-{end}, atom_count={atom_count}")]
    InvalidBondEndpoint {
        bond: BondId,
        begin: AtomId,
        end: AtomId,
        atom_count: usize,
    },
    #[error("bond {bond} stereo {stereo:?} requires two stereo atom references")]
    BondStereoAtomsRequired { bond: BondId, stereo: BondStereo },
    #[error("self-loop bond on bond {bond}: atom {atom}")]
    SelfLoopBond { bond: BondId, atom: AtomId },
    #[error("2D coordinate row count mismatch: rows={rows}, atom_count={atom_count}")]
    CoordinateRowCount { rows: usize, atom_count: usize },
    #[error("3D conformer {conformer} row count mismatch: rows={rows}, atom_count={atom_count}")]
    ConformerRowCount {
        conformer: usize,
        rows: usize,
        atom_count: usize,
    },
    #[error(
        "invalid atom reference in substance group {sgroup:?}: {atom}, atom_count={atom_count}"
    )]
    InvalidSubstanceGroupAtom {
        sgroup: SubstanceGroupId,
        atom: AtomId,
        atom_count: usize,
    },
    #[error(
        "invalid bond reference in substance group {sgroup:?}: {bond}, bond_count={bond_count}"
    )]
    InvalidSubstanceGroupBond {
        sgroup: SubstanceGroupId,
        bond: BondId,
        bond_count: usize,
    },
    #[error("invalid parent reference in substance group {sgroup:?}: {parent:?}")]
    InvalidSubstanceGroupParent {
        sgroup: SubstanceGroupId,
        parent: SubstanceGroupId,
    },
}
