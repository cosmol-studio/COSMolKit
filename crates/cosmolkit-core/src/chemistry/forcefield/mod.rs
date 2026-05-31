//! Source-backed RDKit force-field core primitives.

pub(crate) mod core;
pub(crate) mod crystalff;
pub mod mmff;
pub mod uff;

pub use core::{
    AngleConstraintContrib, AngleConstraintContribs, AngleConstraintContribsParams, DihedralOutput,
    DistanceConstraintContrib, DistanceConstraintContribs, DistanceConstraintContribsParams,
    ForceField, ForceFieldContrib, ForceFieldSnapshot, ForceFieldVec3, PositionConstraintContrib,
    TorsionConstraintContrib, compute_dihedral_from_flat, compute_dihedral_from_points,
    compute_dihedral_from_position_vec, normalize_angle_deg,
};
pub use crystalff::{
    TorsionAngleContribM6, TorsionAngleContribs, TorsionAngleContribsParams, calc_torsion_energy,
    calc_torsion_energy_m6,
};
