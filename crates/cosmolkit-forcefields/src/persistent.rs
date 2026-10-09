//! Owned detached state over the single existing MMFF/UFF force-field kernel.
//! Pinned RDKit: 351f8f378f8ad6bbd517980c38896e66bf907af8.
use crate::kernel::{ForceField, ForceFieldKernelError};
use std::{error::Error, fmt};

/// Stable categories shared by persistent evaluators and their factories.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MolecularForceFieldErrorKind {
    MissingConformer,
    Parameterization,
    InvalidParameterization,
    InvalidAtomIndex,
    InvalidFixedAtom,
    CoordinateCount,
    CoordinateShape,
    NonFiniteCoordinate,
    InvalidTolerance,
    Preparation,
    Construction,
    Initialization,
    Rings,
    Kernel,
}

/// Structured detached error; concrete upstream owner causes stay in the chain.
#[derive(Debug)]
pub struct ForceFieldError {
    kind: MolecularForceFieldErrorKind,
    atom_index: Option<usize>,
    component: Option<usize>,
    actual: Option<usize>,
    expected: Option<usize>,
    source: Option<Box<dyn Error + Send + Sync>>,
}
impl ForceFieldError {
    fn input(
        kind: MolecularForceFieldErrorKind,
        atom_index: Option<usize>,
        component: Option<usize>,
        actual: Option<usize>,
        expected: Option<usize>,
    ) -> Self {
        Self {
            kind,
            atom_index,
            component,
            actual,
            expected,
            source: None,
        }
    }
    fn cause(
        kind: MolecularForceFieldErrorKind,
        source: impl Error + Send + Sync + 'static,
    ) -> Self {
        Self {
            kind,
            atom_index: None,
            component: None,
            actual: None,
            expected: None,
            source: Some(Box::new(source)),
        }
    }
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.kind
    }
    pub fn requested(&self) -> Option<usize> {
        (self.kind == MolecularForceFieldErrorKind::MissingConformer)
            .then_some(self.actual)
            .flatten()
    }
    pub fn atom_index(&self) -> Option<usize> {
        self.atom_index
    }
    pub fn component(&self) -> Option<usize> {
        self.component
    }
    pub fn actual(&self) -> Option<usize> {
        self.actual
    }
    pub fn expected(&self) -> Option<usize> {
        self.expected
    }
}
impl fmt::Display for ForceFieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(source) = &self.source {
            return source.fmt(f);
        }
        write!(
            f,
            "force field {:?}: atom={:?}, component={:?}, actual={:?}, expected={:?}",
            self.kind, self.atom_index, self.component, self.actual, self.expected
        )
    }
}
impl Error for ForceFieldError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn Error + 'static))
    }
}

/// Fixed parameterization and independently owned coordinates, with no live molecule.
pub struct PreparedForceField {
    positions: cosmolkit_model::Conformer3D,
    field: Option<ForceField<'static>>,
}
impl fmt::Debug for PreparedForceField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PreparedForceField")
            .field("positions", &self.positions)
            .finish_non_exhaustive()
    }
}
impl PreparedForceField {
    fn from_borrowed(field: ForceField<'static>, positions: cosmolkit_model::Conformer3D) -> Self {
        // RDKit✔️✔️: RDGeom::PointPtrVect &positions() { return d_positions; }
        // Native ownership: the source field's borrowed position handles are
        // released before retaining independent rows; terms and caches move,
        // never copy/rebuild. No lifetime erasure or self-referential storage.
        Self {
            positions,
            field: Some(field.release_position_borrows()),
        }
    }
    fn with_field<T>(
        &mut self,
        action: impl FnOnce(&mut ForceField<'_>) -> Result<T, ForceFieldKernelError>,
    ) -> Result<T, ForceFieldError> {
        // RDKit✔️✔️: ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
        // FFConvenience.h source rebinding; Native owned rows supply the same
        // atom order. Borrowing costs O(N); the existing release helper retains
        // pointer-vector capacity on this Rust build. Contributions and packed
        // distance matrix retain their allocations, even when action fails.
        let field = self
            .field
            .take()
            .expect("persistent field is restored after every action");
        let mut field = field.clear_positions_and_shorten_lifetime();
        field.positions_mut().extend(
            self.positions
                .coordinates_mut()
                .iter_mut()
                .map(|row| row.as_mut_slice()),
        );
        let result = action(&mut field);
        self.field = Some(field.release_position_borrows());
        result.map_err(|source| {
            let invalid_tolerance = source == ForceFieldKernelError::OptimizerBadTolerance;
            let mut error = ForceFieldError::cause(
                if invalid_tolerance {
                    MolecularForceFieldErrorKind::InvalidTolerance
                } else {
                    MolecularForceFieldErrorKind::Kernel
                },
                source,
            );
            if invalid_tolerance {
                error.component = Some(0);
            }
            error
        })
    }
}

fn selected_index(
    coordinates: &cosmolkit_model::CoordinateBlock,
    id: Option<usize>,
) -> Result<usize, ForceFieldError> {
    // RDKit✔️✔️:   if (d_confs.size() == 0) {
    // RDKit✔️✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto conf : d_confs) {
    // RDKit✔️✔️:     if (conf->getId() == cid) {
    // RDKit✔️✔️:       return *conf;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // ROMol::getConformer: Native approved selection is first stored 3D.
    // O(C) explicit ID lookup, O(1) default, no embedding or AddHs.
    let rows = &coordinates.conformers_3d;
    match id {
        None if !rows.is_empty() => Ok(0),
        Some(id) => rows.iter().position(|row| row.id() == id).ok_or_else(|| {
            ForceFieldError::input(
                MolecularForceFieldErrorKind::MissingConformer,
                None,
                None,
                Some(id),
                None,
            )
        }),
        None => Err(ForceFieldError::input(
            MolecularForceFieldErrorKind::MissingConformer,
            None,
            None,
            None,
            None,
        )),
    }
}
fn validate_rows(rows: &[[f64; 3]], expected: usize) -> Result<(), ForceFieldError> {
    // Native API validation precedes coordinate/cache installation. O(N),
    // fixed three-component Rust rows; Python validates arbitrary shape too.
    if rows.len() != expected {
        return Err(ForceFieldError::input(
            MolecularForceFieldErrorKind::CoordinateCount,
            None,
            None,
            Some(rows.len()),
            Some(expected),
        ));
    }
    for (atom, row) in rows.iter().enumerate() {
        for (component, value) in row.iter().enumerate() {
            if !value.is_finite() {
                return Err(ForceFieldError::input(
                    MolecularForceFieldErrorKind::NonFiniteCoordinate,
                    Some(atom),
                    Some(component),
                    None,
                    None,
                ));
            }
        }
    }
    Ok(())
}

/// Build once from detached MMFF input; topology preparation never mutates a live molecule.
pub fn prepare_mmff_force_field(
    topology: &cosmolkit_model::TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    properties: &cosmolkit_model::MoleculeProperties,
    params: &crate::MmffEvaluationParams,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<PreparedForceField, ForceFieldError> {
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    // RDKit✔️✔️:   res->initialize();
    // MMFF::constructForceField's complete source dispatcher and each helper
    // remain in mmff::builder::construct_force_field_with_props, invoked once.
    // Same MMFF preparation/term ordering and source variant parser; O(N)
    // selected coordinate ownership copy is deliberate Native snapshot cost.
    let props = crate::mmff::mol_properties::MmffMolProperties::new_prepared(
        topology,
        properties.prop("_MMFFSanitized").is_some(),
        &params.mmff_variant,
        0,
        supplied_rings,
    )
    .map_err(|source| {
        ForceFieldError::cause(MolecularForceFieldErrorKind::Parameterization, source)
    })?;
    if !props.is_valid() {
        return Err(ForceFieldError::input(
            MolecularForceFieldErrorKind::InvalidParameterization,
            None,
            None,
            None,
            None,
        ));
    }
    let selected = selected_index(coordinates, params.conformer_id)?;
    let stored = &coordinates.conformers_3d;
    let row = &stored[selected];
    validate_rows(row.coordinates(), topology.atoms.len())?;
    let mut positions = row.coordinates().to_vec();
    let context = crate::mmff::builder::MmffConformerContext {
        two_d: &coordinates.conformers_2d,
        before: &stored[..selected],
        selected_id: row.id(),
        selected_is_3d: row.is_3d(),
        selected_props: row.props(),
        after: &stored[selected + 1..],
        source_dimension: coordinates.source_coordinate_dim,
        source_order: coordinates.source_conformer_order.as_deref(),
    };
    let mut field = crate::mmff::builder::construct_force_field_with_props(
        &props.topology,
        &props,
        positions.iter_mut().map(|row| row.as_mut_slice()).collect(),
        &context,
        properties,
        params.non_bonded_threshold,
        params.ignore_interfragment_interactions,
    )
    .map_err(|source| ForceFieldError::cause(MolecularForceFieldErrorKind::Construction, source))?;
    field.initialize().map_err(|source| {
        ForceFieldError::cause(MolecularForceFieldErrorKind::Initialization, source)
    })?;
    let field = field.release_position_borrows();
    Ok(PreparedForceField::from_borrowed(
        field,
        cosmolkit_model::Conformer3D::new(row.id(), positions, row.is_3d()),
    ))
}

/// Build UFF once from detached prepared valence/ring state and an existing 3D row.
#[allow(clippy::too_many_arguments)]
pub fn prepare_uff_force_field(
    topology: &cosmolkit_model::TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    valence: &cosmolkit_core::ValenceAssignment,
    rings: &cosmolkit_core::RingInfo,
    properties: &cosmolkit_model::MoleculeProperties,
    params: &crate::UffEvaluationParams,
) -> Result<PreparedForceField, ForceFieldError> {
    // RDKit✔️✔️:   bool foundAll;
    // RDKit✔️✔️:   AtomicParamVect params;
    // RDKit✔️✔️:   boost::tie(params, foundAll) = getAtomTypes(mol);
    // RDKit✔️✔️:   return constructForceField(mol, params, vdwThresh, confId,
    // RDKit✔️✔️:                              ignoreInterfragInteractions);
    // Complete source dispatch and getAtomTypes helper remain inline in the
    // unique uff::builder owners. Native owned selection clones only one row;
    // O(N) ownership cost once, no topology or unrelated conformer copy.
    let prepared =
        crate::uff::api::prepare_parameter_query(topology, valence).map_err(|source| {
            ForceFieldError::cause(MolecularForceFieldErrorKind::Preparation, source)
        })?;
    let index = selected_index(coordinates, params.conformer_id)?;
    let rows = &coordinates.conformers_3d;
    let row = &rows[index];
    validate_rows(row.coordinates(), topology.atoms.len())?;
    let mut selected = row.clone();
    let context = crate::uff::builder::UffConformerContext {
        two_d: &coordinates.conformers_2d,
        before: &rows[..index],
        selected_id: row.id(),
        selected_is_3d: row.is_3d(),
        selected_props: row.props(),
        after: &rows[index + 1..],
        source_dimension: coordinates.source_coordinate_dim,
        source_order: coordinates.source_conformer_order.as_deref(),
    };
    let mut diagnostics = Vec::new();
    let mut field = crate::uff::builder::construct_force_field_with_automatic_typing_from_selected(
        topology,
        &mut selected,
        &context,
        prepared.typing_state,
        rings,
        valence,
        properties,
        &mut diagnostics,
        params.vdw_threshold,
        params.ignore_interfragment_interactions,
    )
    .map_err(|source| ForceFieldError::cause(MolecularForceFieldErrorKind::Construction, source))?;
    field.initialize().map_err(|source| {
        ForceFieldError::cause(MolecularForceFieldErrorKind::Initialization, source)
    })?;
    let field = field.release_position_borrows();
    Ok(PreparedForceField::from_borrowed(field, selected))
}

impl PreparedForceField {
    pub fn position(&self, atom: usize) -> Result<[f64; 3], ForceFieldError> {
        self.positions
            .coordinates()
            .get(atom)
            .copied()
            .ok_or_else(|| {
                ForceFieldError::input(
                    MolecularForceFieldErrorKind::InvalidAtomIndex,
                    Some(atom),
                    None,
                    Some(atom),
                    Some(self.positions.coordinates().len()),
                )
            })
    }
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.positions.coordinates().to_vec()
    }
    pub fn set_position(&mut self, atom: usize, row: [f64; 3]) -> Result<(), ForceFieldError> {
        // Native API: validate every input before touching state; fixed atoms
        // may be dragged to a new anchor. O(1) row write, no coordinate clone.
        self.position(atom)?;
        validate_rows(std::slice::from_ref(&row), 1).map_err(|mut error| {
            error.atom_index = Some(atom);
            error
        })?;
        self.field
            .as_mut()
            .expect("field available")
            .init_distance_matrix()
            .map_err(|source| {
                ForceFieldError::cause(MolecularForceFieldErrorKind::Kernel, source)
            })?;
        self.positions.coordinates_mut()[atom] = row;
        Ok(())
    }
    pub fn set_positions(&mut self, rows: &[[f64; 3]]) -> Result<(), ForceFieldError> {
        // Native API: one O(N) validation pass then O(N) copy into existing
        // storage. Only the source distance matrix is invalidated.
        validate_rows(rows, self.positions.coordinates().len())?;
        self.field
            .as_mut()
            .expect("field available")
            .init_distance_matrix()
            .map_err(|source| {
                ForceFieldError::cause(MolecularForceFieldErrorKind::Kernel, source)
            })?;
        self.positions.coordinates_mut().copy_from_slice(rows);
        Ok(())
    }
    pub fn fixed_atoms(&self) -> Vec<usize> {
        self.field
            .as_ref()
            .expect("field available")
            .fixed_points()
            .iter()
            .map(|&atom| atom as usize)
            .collect()
    }
    pub fn set_fixed_atoms(&mut self, atoms: &[usize]) -> Result<(), ForceFieldError> {
        // RDKit✔️✔️: INT_VECT &fixedPoints() { return d_fixedPoints; }
        // Native complete-set replacement, atomic validation and canonical
        // sorted/deduplicated set. Fixed points still contribute to energies.
        for &atom in atoms {
            self.position(atom).map_err(|mut error| {
                error.kind = MolecularForceFieldErrorKind::InvalidFixedAtom;
                error
            })?;
        }
        let mut fixed = atoms
            .iter()
            .map(|&atom| {
                i32::try_from(atom).map_err(|_| {
                    ForceFieldError::input(
                        MolecularForceFieldErrorKind::InvalidFixedAtom,
                        Some(atom),
                        None,
                        Some(atom),
                        Some(i32::MAX as usize),
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        fixed.sort_unstable();
        fixed.dedup();
        *self
            .field
            .as_mut()
            .expect("field available")
            .fixed_points_mut() = fixed;
        Ok(())
    }
}

impl PreparedForceField {
    pub fn energy(&mut self) -> Result<f64, ForceFieldError> {
        // RDKit✔️✔️:     double e = d_contrib->getEnergy(pos);
        // RDKit✔️✔️:     res += e;
        // Exact complete source body remains in calc_energy_current, including
        // scatter and original contribution order. No fixed-point energy mask.
        self.with_field(|field| field.calc_energy_current(None))
    }
    pub fn gradient(&mut self) -> Result<Vec<[f64; 3]>, ForceFieldError> {
        // RDKit✔️✔️:     d_contrib->getGrad(pos, grad);
        // RDKit✔️✔️:       grad[idx + di] = 0.0;
        // Delegate accumulation/fixed zeroing to the unique calc_grad_current.
        // O(N) independent result, no mutable view of owned live coordinates.
        self.gradient_with_fixed_mask(true)
    }

    /// Full energy derivative, including fixed-atom rows. Physical force is
    /// its negative. Only the fixed-point mask is omitted; energy terms and
    /// the pinned set are unchanged.
    pub fn gradient_unconstrained(&mut self) -> Result<Vec<[f64; 3]>, ForceFieldError> {
        self.gradient_with_fixed_mask(false)
    }

    fn gradient_with_fixed_mask(
        &mut self,
        apply_fixed_mask: bool,
    ) -> Result<Vec<[f64; 3]>, ForceFieldError> {
        let count = self.positions.coordinates().len();
        self.with_field(|field| {
            let mut gradient = vec![0.0; 3 * count];
            if apply_fixed_mask {
                field.calc_grad_current(&mut gradient)?;
            } else {
                field.calc_grad_current_unconstrained(&mut gradient)?;
            }
            Ok(gradient
                .chunks_exact(3)
                .map(|row| [row[0], row[1], row[2]])
                .collect())
        })
    }
    pub fn energy_gradient(&mut self) -> Result<(f64, Vec<[f64; 3]>), ForceFieldError> {
        // RDKit✔️✔️:     double e = d_contrib->getEnergy(pos);
        // RDKit✔️✔️:     d_contrib->getGrad(pos, grad);
        // Both source owners share the same field/cache and current rows.
        // One binding pass; original energy and gradient accumulation are
        // separate source calls. O(N) result, no copying/rebuilding terms.
        let count = self.positions.coordinates().len();
        self.with_field(|field| {
            let energy = field.calc_energy_current(None)?;
            let mut gradient = vec![0.0; 3 * count];
            field.calc_grad_current(&mut gradient)?;
            Ok((
                energy,
                gradient
                    .chunks_exact(3)
                    .map(|row| [row[0], row[1], row[2]])
                    .collect(),
            ))
        })
    }
}

impl PreparedForceField {
    pub fn minimize(
        &mut self,
        max_iterations: u32,
        force_tolerance: f64,
        energy_tolerance: f64,
    ) -> Result<(bool, u32, f64), ForceFieldError> {
        // RDKit✔️✔️:   unsigned int numIters = 0;
        // RDKit✔️✔️:       BFGSOpt::minimize(dim, points.data(), forceTol, numIters, finalForce, eCalc,
        // RDKit✔️✔️:                         gCalc, snapshotFreq, snapshotVect, energyTol, maxIts);
        // RDKit✔️✔️:   this->gather(points.data());
        // Same dense BFGS owner starts fresh history each call. Native outcome
        // exposes its existing counter. Forward energyTol unchanged: the pinned
        // BFGS implementation explicitly ignores funcTol (RDUNUSED_PARAM).
        self.with_field(|field| {
            let outcome =
                field.minimize_with_iterations(max_iterations, force_tolerance, energy_tolerance);
            // A failed callback can cache trial coordinates without gathering;
            // invalidate those distances even on Err before the next query.
            field.init_distance_matrix()?;
            let (status, iterations) = outcome?;
            let energy = field.calc_energy_current(None)?;
            Ok((status == 0, iterations, energy))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer3D, CoordinateBlock,
        Element, Hybridization, MoleculeProperties, TopologyBlock,
    };
    fn pair(bonded: bool, distance: f64) -> (TopologyBlock, CoordinateBlock) {
        let atoms = (0..2)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                )
            })
            .collect();
        let bonds = if bonded {
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )]
        } else {
            vec![]
        };
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let mut coordinates = CoordinateBlock::default();
        coordinates.conformers_3d.push(Conformer3D::new(
            7,
            vec![[0., 0., 0.], [distance, 0., 0.]],
            true,
        ));
        (topology, coordinates)
    }
    fn mmff(topology: &TopologyBlock, coordinates: &CoordinateBlock) -> PreparedForceField {
        prepare_mmff_force_field(
            topology,
            coordinates,
            &MoleculeProperties::default(),
            &crate::MmffEvaluationParams::default(),
            None,
        )
        .unwrap()
    }
    fn uff(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        threshold: f64,
        ignore: bool,
    ) -> PreparedForceField {
        let valence = cosmolkit_core::assign_valence_for_topology(
            topology,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .unwrap();
        let rings = cosmolkit_core::symmetrized_sssr(
            topology,
            &cosmolkit_core::RingSearchParams::default(),
        )
        .unwrap();
        prepare_uff_force_field(
            topology,
            coordinates,
            &valence,
            &rings,
            &MoleculeProperties::default(),
            &crate::UffEvaluationParams {
                vdw_threshold: threshold,
                ignore_interfragment_interactions: ignore,
                ..Default::default()
            },
        )
        .unwrap()
    }
    #[test]
    fn persistent_mmff_evaluation_and_minimize_reuse_original_owners() {
        let (topology, coordinates) = pair(true, 2.5);
        let expected = crate::evaluate_mmff(
            &topology,
            &coordinates,
            &MoleculeProperties::default(),
            &crate::MmffEvaluationParams::default(),
            None,
        )
        .unwrap()
        .unwrap();
        let mut field = mmff(&topology, &coordinates);
        let (energy, gradient) = field.energy_gradient().unwrap();
        assert_eq!(energy, expected.energy);
        assert_eq!(
            gradient.iter().flatten().copied().collect::<Vec<_>>(),
            expected.gradient
        );
        assert_eq!(field.gradient().unwrap(), gradient);
        let mut work_topology = topology.clone();
        let mut work_coordinates = coordinates.clone();
        let expected = crate::optimize_mmff_single(
            &mut work_topology,
            &mut work_coordinates,
            &MoleculeProperties::default(),
            &crate::MmffOptimizationParams {
                max_iterations: 20,
                ..Default::default()
            },
            None,
        )
        .unwrap();
        let (converged, iterations, energy) = field.minimize(20, 1.0e-4, 1.0e-6).unwrap();
        assert_eq!(converged, expected.status == 0);
        assert!(iterations > 0 && iterations <= 20);
        assert_eq!(
            field.positions(),
            work_coordinates.conformers_3d[0].coordinates()
        );
        assert_eq!(energy, field.energy().unwrap());
        assert_eq!(coordinates.conformers_3d[0].coordinates()[1], [2.5, 0., 0.]);
    }
    #[test]
    fn persistent_uff_cache_refresh_and_independent_snapshots() {
        let (topology, coordinates) = pair(true, 2.5);
        let mut field = uff(&topology, &coordinates, 10., true);
        let initial = field.energy().unwrap();
        let mut snapshot = field.positions();
        snapshot[1][0] = 1.7;
        assert_eq!(field.energy().unwrap(), initial);
        field.set_position(1, snapshot[1]).unwrap();
        let changed = field.energy().unwrap();
        assert_ne!(initial, changed);
        let mut current = coordinates.clone();
        current.conformers_3d[0]
            .coordinates_mut()
            .copy_from_slice(&snapshot);
        assert_eq!(
            changed,
            uff(&topology, &current, 10., true).energy().unwrap()
        );
        assert_eq!(
            field.gradient().unwrap(),
            uff(&topology, &current, 10., true).gradient().unwrap()
        );
        assert_eq!(coordinates.conformers_3d[0].coordinates()[1], [2.5, 0., 0.]);
    }
    #[test]
    fn persistent_invalid_updates_are_atomic_and_errors_structured() {
        let (topology, coordinates) = pair(true, 2.5);
        let mut field = mmff(&topology, &coordinates);
        field.set_fixed_atoms(&[0]).unwrap();
        let before = field.positions();
        let energy = field.energy().unwrap();
        assert_eq!(
            field.set_position(2, [0.; 3]).unwrap_err().kind(),
            MolecularForceFieldErrorKind::InvalidAtomIndex
        );
        let error = field.set_position(1, [1., f64::NAN, 3.]).unwrap_err();
        assert_eq!(
            (error.kind(), error.atom_index(), error.component()),
            (
                MolecularForceFieldErrorKind::NonFiniteCoordinate,
                Some(1),
                Some(1)
            )
        );
        assert_eq!(
            field.set_positions(&[[0.; 3]]).unwrap_err().kind(),
            MolecularForceFieldErrorKind::CoordinateCount
        );
        assert_eq!(
            field
                .set_positions(&[[0.; 3], [f64::INFINITY, 0., 0.]])
                .unwrap_err()
                .kind(),
            MolecularForceFieldErrorKind::NonFiniteCoordinate
        );
        assert_eq!(
            field.set_fixed_atoms(&[1, 2]).unwrap_err().kind(),
            MolecularForceFieldErrorKind::InvalidFixedAtom
        );
        assert_eq!(field.positions(), before);
        assert_eq!(field.fixed_atoms(), [0]);
        assert_eq!(field.energy().unwrap(), energy);
    }
    #[test]
    fn persistent_fixed_set_replacement_and_explicit_anchor_update() {
        let (topology, coordinates) = pair(true, 2.5);
        let mut field = mmff(&topology, &coordinates);
        field.set_fixed_atoms(&[0, 0]).unwrap();
        assert_eq!(field.fixed_atoms(), [0]);
        field.set_position(0, [0.2, 0.1, 0.]).unwrap();
        let anchor = field.position(0).unwrap();
        assert_eq!(field.gradient().unwrap()[0], [0.; 3]);
        assert_ne!(field.gradient().unwrap()[1], [0.; 3]);
        field.minimize(20, 1.0e-4, 1.0e-6).unwrap();
        assert_eq!(field.position(0).unwrap(), anchor);
        field.set_fixed_atoms(&[1]).unwrap();
        assert_eq!(field.fixed_atoms(), [1]);
        field.set_fixed_atoms(&[]).unwrap();
        assert!(field.fixed_atoms().is_empty());
    }
    #[test]
    fn persistent_unconstrained_gradient_preserves_pins_and_exact_accumulation() {
        let (topology, coordinates) = pair(true, 2.5);
        for is_mmff in [false, true] {
            let make = || {
                if is_mmff {
                    mmff(&topology, &coordinates)
                } else {
                    uff(&topology, &coordinates, 10., true)
                }
            };
            let mut pinned = make();
            let mut free = make();
            pinned.set_fixed_atoms(&[0]).unwrap();
            for x in [2.5, 1.9] {
                pinned.set_position(1, [x, 0.2, 0.]).unwrap();
                free.set_position(1, [x, 0.2, 0.]).unwrap();
                let before = pinned.positions();
                let expected = free.gradient().unwrap();
                let actual = pinned.gradient_unconstrained().unwrap();
                assert_eq!(
                    actual
                        .iter()
                        .flatten()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>(),
                    expected
                        .iter()
                        .flatten()
                        .map(|x| x.to_bits())
                        .collect::<Vec<_>>()
                );
                assert_ne!(actual[0], [0.; 3]);
                let constrained = pinned.gradient().unwrap();
                assert_eq!(constrained[0], [0.; 3]);
                assert_eq!(
                    constrained[1].map(f64::to_bits),
                    actual[1].map(f64::to_bits)
                );
                assert_eq!(
                    pinned.energy().unwrap().to_bits(),
                    free.energy().unwrap().to_bits()
                );
                assert_eq!(pinned.positions(), before);
                assert_eq!(pinned.fixed_atoms(), [0]);
            }
            let anchor = pinned.position(0).unwrap();
            pinned.minimize(2, 1e-4, 1e-6).unwrap();
            assert_eq!(pinned.position(0).unwrap(), anchor);
        }
    }
    #[test]
    fn persistent_repeated_minimize_restarts_history_at_current_coordinates() {
        let (topology, coordinates) = pair(true, 2.5);
        let mut field = mmff(&topology, &coordinates);
        assert_eq!(field.minimize(0, 1.0e-4, 1.0e-6).unwrap().0, false);
        assert_eq!(field.minimize(0, 1.0e-4, 1.0e-6).unwrap().1, 0);
        field.minimize(2, 1.0e-4, 1.0e-6).unwrap();
        let mut current = coordinates.clone();
        current.conformers_3d[0]
            .coordinates_mut()
            .copy_from_slice(&field.positions());
        let mut fresh = mmff(&topology, &current);
        let repeated = field.minimize(2, 1.0e-4, 1.0e-6).unwrap();
        assert_eq!(repeated, fresh.minimize(2, 1.0e-4, 1.0e-6).unwrap());
        assert_eq!(field.positions(), fresh.positions());
        let before = field.positions();
        let error = field.minimize(20, 0., 1.0e-6).unwrap_err();
        assert_eq!(error.kind(), MolecularForceFieldErrorKind::InvalidTolerance);
        assert!(error.source().is_some());
        assert_eq!(field.positions(), before);
        assert_eq!(field.energy().unwrap(), fresh.energy().unwrap());
    }
    #[test]
    fn persistent_energy_tolerance_preserves_pinned_bfgs_semantics() {
        // RDKit BFGSOpt.h:184-189:
        // RDKit✔️✔️:   RDUNUSED_PARAM(funcTol);
        // The public parameter is forwarded, not a CK-only stopping rule.
        let (topology, coordinates) = pair(true, 2.5);
        let mut baseline = mmff(&topology, &coordinates);
        let expected = baseline.minimize(2, 1e-4, 1e-6).unwrap();
        for energy_tolerance in [0.0, 1e-12, 1e6] {
            let mut field = mmff(&topology, &coordinates);
            let actual = field.minimize(2, 1e-4, energy_tolerance).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(field.positions(), baseline.positions());
            assert_eq!(field.gradient().unwrap(), baseline.gradient().unwrap());
        }
    }
    #[test]
    fn persistent_nonbonded_selection_is_frozen_at_construction() {
        let (topology, coordinates) = pair(false, 100.);
        let mut field = uff(&topology, &coordinates, 1., false);
        assert_eq!(field.energy().unwrap(), 0.);
        field.set_position(1, [2., 0., 0.]).unwrap();
        assert_eq!(field.energy().unwrap(), 0.);
        let mut current = coordinates.clone();
        current.conformers_3d[0].coordinates_mut()[1] = [2., 0., 0.];
        assert_ne!(uff(&topology, &current, 1., false).energy().unwrap(), 0.);
        assert_eq!(field.minimize(20, 1.0e-4, 1.0e-6).unwrap(), (true, 0, 0.));
    }
    #[test]
    fn persistent_owned_handle_outlives_inputs_and_selects_existing_3d_id() {
        let mut field = {
            let (topology, mut coordinates) = pair(true, 2.5);
            coordinates.conformers_3d.push(Conformer3D::new(
                42,
                vec![[0., 0., 0.], [1.8, 0., 0.]],
                true,
            ));
            let selected = prepare_mmff_force_field(
                &topology,
                &coordinates,
                &MoleculeProperties::default(),
                &crate::MmffEvaluationParams {
                    conformer_id: Some(42),
                    ..Default::default()
                },
                None,
            )
            .unwrap();
            assert_eq!(selected.position(1).unwrap(), [1.8, 0., 0.]);
            mmff(&topology, &coordinates)
        };
        assert_eq!(field.position(1).unwrap(), [2.5, 0., 0.]);
        assert!(field.energy().unwrap().is_finite());
        field.set_positions(&[[0., 0., 0.], [2.1, 0., 0.]]).unwrap();
        assert_eq!(field.position(1).unwrap(), [2.1, 0., 0.]);
    }
    #[test]
    fn persistent_mmff_invalid_parameterization_fails_construction() {
        let topology = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::HE))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let mut coordinates = CoordinateBlock::default();
        coordinates
            .conformers_3d
            .push(Conformer3D::new(0, vec![[0.; 3]], true));
        let error = prepare_mmff_force_field(
            &topology,
            &coordinates,
            &MoleculeProperties::default(),
            &crate::MmffEvaluationParams::default(),
            None,
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            MolecularForceFieldErrorKind::InvalidParameterization
        );
        let (topology, _) = pair(true, 2.5);
        let error = prepare_mmff_force_field(
            &topology,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
            &crate::MmffEvaluationParams::default(),
            None,
        )
        .unwrap_err();
        assert_eq!(error.kind(), MolecularForceFieldErrorKind::MissingConformer);
    }
}
