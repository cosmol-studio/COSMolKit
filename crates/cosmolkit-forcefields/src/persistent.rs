//! Owned detached state over the single existing MMFF/UFF force-field kernel.
//! Pinned RDKit: 351f8f378f8ad6bbd517980c38896e66bf907af8.
use crate::kernel::{
    AngleConstraintContrib, DistanceConstraintContrib, ForceField, ForceFieldKernelError,
    PositionConstraintContrib, TorsionConstraintContrib,
};
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
    InvalidConstraintBounds,
    InvalidConstraintRange,
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
#[derive(Clone, Copy)]
enum ConstraintIndexKind {
    Distance,
    Angle,
    Torsion,
}
impl fmt::Debug for PreparedForceField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PreparedForceField")
            .field("positions", &self.positions)
            .finish_non_exhaustive()
    }
}
impl PreparedForceField {
    pub fn add_position_constraint(
        &mut self,
        atom: usize,
        max_displacement: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldError> {
        // RDKit✔️✔️:   auto *constraint = new UFF::PositionConstraintContrib(
        // RDKit✔️✔️:       self->field.get(), idx, maxDispl, forceConstant);
        // RDKit✔️✔️:   self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
        // One current Point3 copy; no rebind, numeric filter or coordinate write.
        let constraint = PositionConstraintContrib::new(
            self.positions.coordinates(),
            self.constraint_index(atom, 0, ConstraintIndexKind::Distance)?,
            max_displacement,
            force_constant,
        )
        .map_err(constraint_error)?;
        self.field
            .as_mut()
            .expect("field available")
            .add_contribution(Box::new(constraint));
        Ok(())
    }
    pub fn add_torsion_constraint(
        &mut self,
        atom1: usize,
        atom2: usize,
        atom3: usize,
        atom4: usize,
        relative: bool,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldError> {
        // RDKit✔️✔️:   auto *constraint = new UFF::TorsionConstraintContrib(
        // RDKit✔️✔️:       self->field.get(), idx1, idx2, idx3, idx4, relative, minDihedralDeg,
        // RDKit✔️✔️:       maxDihedralDeg, forceConstant);
        // RDKit✔️✔️:   self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
        // The one constructor captures only these four current rows and normalizes
        // source bounds before append. No cache reset, term copies or rebinding.
        if min_dihedral_deg > max_dihedral_deg {
            return Err(constraint_error(ForceFieldKernelError::TorsionBoundsOrder));
        }
        let constraint = TorsionConstraintContrib::new_relative(
            self.positions.coordinates(),
            self.constraint_index(atom1, 0, ConstraintIndexKind::Torsion)?,
            self.constraint_index(atom2, 1, ConstraintIndexKind::Torsion)?,
            self.constraint_index(atom3, 2, ConstraintIndexKind::Torsion)?,
            self.constraint_index(atom4, 3, ConstraintIndexKind::Torsion)?,
            relative,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        )
        .map_err(constraint_error)?;
        self.field
            .as_mut()
            .expect("field available")
            .add_contribution(Box::new(constraint));
        Ok(())
    }
    pub fn add_angle_constraint(
        &mut self,
        atom1: usize,
        atom2: usize,
        atom3: usize,
        relative: bool,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldError> {
        // RDKit✔️✔️:   auto *constraint = new UFF::AngleConstraintContrib(
        // RDKit✔️✔️:       self->field.get(), idx1, idx2, idx3, relative, minAngleDeg, maxAngleDeg,
        // RDKit✔️✔️:       forceConstant);
        // RDKit✔️✔️:   self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
        // The shared constructor reads only these three current owned rows.
        let constraint = AngleConstraintContrib::new_relative(
            self.positions.coordinates(),
            self.constraint_index(atom1, 0, ConstraintIndexKind::Angle)?,
            self.constraint_index(atom2, 1, ConstraintIndexKind::Angle)?,
            self.constraint_index(atom3, 2, ConstraintIndexKind::Angle)?,
            relative,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        )
        .map_err(constraint_error)?;
        self.field
            .as_mut()
            .expect("field available")
            .add_contribution(Box::new(constraint));
        Ok(())
    }
    pub fn add_distance_constraint(
        &mut self,
        atom1: usize,
        atom2: usize,
        relative: bool,
        min_distance: f64,
        max_distance: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldError> {
        // RDKit✔️✔️:   auto *constraint = new UFF::DistanceConstraintContrib(
        // RDKit✔️✔️:       self->field.get(), idx1, idx2, relative, minLen, maxLen, forceConstant);
        // RDKit✔️✔️:   self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
        // MMFF aliases the same contribution. Indexed owned rows give constant
        // construction work; no N-row rebinding, clone, cache reset or rebuild.
        let constraint = DistanceConstraintContrib::new_relative(
            self.positions.coordinates(),
            self.constraint_index(atom1, 0, ConstraintIndexKind::Distance)?,
            self.constraint_index(atom2, 1, ConstraintIndexKind::Distance)?,
            relative,
            min_distance,
            max_distance,
            force_constant,
        )
        .map_err(constraint_error)?;
        self.field
            .as_mut()
            .expect("field available")
            .add_contribution(Box::new(constraint));
        Ok(())
    }

    fn constraint_index(
        &self,
        atom: usize,
        component: usize,
        kind: ConstraintIndexKind,
    ) -> Result<u32, ForceFieldError> {
        use crate::kernel::{AngleIndexArgument, ForceFieldIndexArgument, TorsionIndexArgument};
        // Native identifier representation boundary: do not narrow a usize ID
        // into a different valid source u32 index. Constant work, no mutation.
        // Validate each source range in argument order before converting the
        // next ID, retaining the concrete original cause for representable IDs.
        let count = self.positions.coordinates().len();
        if atom >= count {
            if let Ok(index) = u32::try_from(atom) {
                let source = match kind {
                    ConstraintIndexKind::Distance => ForceFieldKernelError::IndexOutOfRange {
                        argument: if component == 0 {
                            ForceFieldIndexArgument::I
                        } else {
                            ForceFieldIndexArgument::J
                        },
                        index,
                        upper_bound: count as u32,
                    },
                    ConstraintIndexKind::Angle => ForceFieldKernelError::AngleIndexOutOfRange {
                        argument: [
                            AngleIndexArgument::First,
                            AngleIndexArgument::Second,
                            AngleIndexArgument::Third,
                        ][component],
                        index,
                        upper_bound: count,
                    },
                    ConstraintIndexKind::Torsion => ForceFieldKernelError::TorsionIndexOutOfRange {
                        argument: [
                            TorsionIndexArgument::First,
                            TorsionIndexArgument::Second,
                            TorsionIndexArgument::Third,
                            TorsionIndexArgument::Fourth,
                        ][component],
                        index,
                        upper_bound: count,
                    },
                };
                return Err(constraint_error(source));
            }
        }
        u32::try_from(atom).map_err(|_| {
            ForceFieldError::input(
                MolecularForceFieldErrorKind::InvalidAtomIndex,
                Some(atom),
                Some(component),
                Some(atom),
                Some(self.positions.coordinates().len()),
            )
        })
    }
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

fn constraint_error(source: ForceFieldKernelError) -> ForceFieldError {
    use crate::kernel::ForceFieldIndexArgument;
    let mut error = ForceFieldError::cause(MolecularForceFieldErrorKind::Kernel, source);
    match source {
        ForceFieldKernelError::IndexOutOfRange {
            argument,
            index,
            upper_bound,
        } => {
            error.kind = MolecularForceFieldErrorKind::InvalidAtomIndex;
            error.atom_index = Some(index as usize);
            error.actual = Some(index as usize);
            error.expected = Some(upper_bound as usize);
            error.component = Some(match argument {
                ForceFieldIndexArgument::I => 0,
                ForceFieldIndexArgument::J => 1,
            });
        }
        ForceFieldKernelError::BadBounds => {
            error.kind = MolecularForceFieldErrorKind::InvalidConstraintBounds
        }
        ForceFieldKernelError::TorsionBoundsOrder => {
            error.kind = MolecularForceFieldErrorKind::InvalidConstraintBounds
        }
        ForceFieldKernelError::TorsionIndexOutOfRange {
            argument,
            index,
            upper_bound,
        } => {
            error.kind = MolecularForceFieldErrorKind::InvalidAtomIndex;
            error.atom_index = Some(index as usize);
            error.actual = Some(index as usize);
            error.expected = Some(upper_bound);
            error.component = Some(match argument {
                crate::kernel::TorsionIndexArgument::First => 0,
                crate::kernel::TorsionIndexArgument::Second => 1,
                crate::kernel::TorsionIndexArgument::Third => 2,
                crate::kernel::TorsionIndexArgument::Fourth => 3,
            });
        }
        ForceFieldKernelError::AngleBoundsOrder => {
            error.kind = MolecularForceFieldErrorKind::InvalidConstraintBounds
        }
        ForceFieldKernelError::AngleOutOfRange { bound } => {
            error.kind = MolecularForceFieldErrorKind::InvalidConstraintRange;
            error.component = Some(match bound {
                crate::kernel::AngleRangeBound::Minimum => 0,
                crate::kernel::AngleRangeBound::Maximum => 1,
            });
        }
        ForceFieldKernelError::AngleIndexOutOfRange {
            argument,
            index,
            upper_bound,
        } => {
            error.kind = MolecularForceFieldErrorKind::InvalidAtomIndex;
            error.atom_index = Some(index as usize);
            error.actual = Some(index as usize);
            error.expected = Some(upper_bound);
            error.component = Some(match argument {
                crate::kernel::AngleIndexArgument::First => 0,
                crate::kernel::AngleIndexArgument::Second => 1,
                crate::kernel::AngleIndexArgument::Third => 2,
            });
        }
        _ => {}
    }
    error
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
        params.mmff_variant.as_rdkit_str(),
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
    fn constraint_only(rows: Vec<[f64; 3]>) -> PreparedForceField {
        let mut borrowed = rows.clone();
        let mut field = ForceField::new(3);
        field
            .positions_mut()
            .extend(borrowed.iter_mut().map(|row| row.as_mut_slice()));
        field.initialize().unwrap();
        PreparedForceField::from_borrowed(
            field.release_position_borrows(),
            Conformer3D::new(0, rows, true),
        )
    }
    #[test]
    fn persistent_position_constraint_keeps_original_anchor_and_source_numbers() {
        let mut field = constraint_only(vec![[1., 2., 3.]]);
        field.add_position_constraint(0, 1., 4.).unwrap();
        for point in [[1., 2., 3.], [2., 2., 3.]] {
            field.set_position(0, point).unwrap();
            assert_eq!(field.energy().unwrap().to_bits(), 0.0f64.to_bits());
            assert_eq!(field.gradient().unwrap(), [[0.; 3]]);
        }
        field.set_position(0, [3., 2., 3.]).unwrap();
        assert_eq!(field.energy().unwrap().to_bits(), 2.0f64.to_bits());
        assert_eq!(field.gradient().unwrap(), [[4., 0., 0.]]);
        field.add_position_constraint(0, 0., 2.).unwrap();
        field.set_position(0, [1., 2., 3.]).unwrap();
        assert_eq!(field.energy().unwrap().to_bits(), 4.0f64.to_bits());
        assert_eq!(field.gradient().unwrap(), [[-4., 0., 0.]]);
        let before = field.energy().unwrap();
        assert_eq!(
            field.add_position_constraint(8, 0., 2.).unwrap_err().kind(),
            MolecularForceFieldErrorKind::InvalidAtomIndex
        );
        assert_eq!(field.energy().unwrap().to_bits(), before.to_bits());
        // Source permits negative displacement/force and NaN parameters.
        let mut negative = constraint_only(vec![[0.; 3]]);
        negative.add_position_constraint(0, -1., -4.).unwrap();
        assert_eq!(negative.energy().unwrap().to_bits(), (-2.0f64).to_bits());
        assert_eq!(negative.gradient().unwrap(), [[0.; 3]]);
        negative.add_position_constraint(0, f64::NAN, 1.).unwrap();
        assert!(negative.energy().unwrap().is_nan());
    }
    #[test]
    fn persistent_mixed_constraints_append_order_tie_and_fixed_mask() {
        let mut order = constraint_only(vec![[0.; 3], [3., 0., 0.]]);
        for force in [2.0e16, -2.0e16, 2.] {
            order
                .add_distance_constraint(0, 1, false, 2., 2., force)
                .unwrap();
        }
        assert_eq!(order.energy().unwrap().to_bits(), 1.0f64.to_bits());
        order.set_fixed_atoms(&[0]).unwrap();
        assert_eq!(order.gradient().unwrap()[0], [0.; 3]);
        assert_eq!(order.gradient_unconstrained().unwrap()[0], [-2., 0., 0.]);
        let rows = vec![[1., 0., 0.], [0.; 3], [0., 1., 0.], [1., 1., 0.]];
        let mut tie = constraint_only(rows);
        tie.add_torsion_constraint(0, 1, 2, 3, false, 170., 190., 2.)
            .unwrap();
        // Source chooses maximum on equal absolute endpoint distances.
        assert_eq!(tie.energy().unwrap().to_bits(), 57800.0f64.to_bits());
        assert_eq!(tie.gradient().unwrap()[0][2].to_bits(), 0x40e3062429863ecf);
        tie.add_angle_constraint(0, 1, 2, false, f64::NAN, 100., 1.)
            .unwrap();
        tie.add_torsion_constraint(0, 1, 2, 3, false, f64::NAN, 0., 1.)
            .unwrap();
        assert_eq!(tie.energy().unwrap().to_bits(), 57800.0f64.to_bits());
    }
    #[test]
    fn persistent_mixed_constraints_three_models_and_active_minimization() {
        let (topology, coordinates) = pair(true, 2.5);
        for (variant, expected_energy, expected_gradient, final_energy) in [
            (
                "UFF",
                0x407c491f6ef37e79,
                0x408a3e614ec623d7,
                0x402b8181218a4b74,
            ),
            (
                "MMFF94",
                0x407fc2d63b90492d,
                0x409b326feceb2652,
                0x402af985a2348490,
            ),
            (
                "MMFF94s",
                0x407fc2d63b90492d,
                0x409b326feceb2652,
                0x402af985a2348490,
            ),
        ] {
            let mut field = if variant == "UFF" {
                uff(&topology, &coordinates, 10., true)
            } else {
                prepare_mmff_force_field(
                    &topology,
                    &coordinates,
                    &MoleculeProperties::default(),
                    &crate::MmffEvaluationParams {
                        mmff_variant: crate::MmffVariant::from_rdkit_constructor_arg(variant),
                        ..Default::default()
                    },
                    None,
                )
                .unwrap()
            };
            field
                .add_distance_constraint(0, 1, false, 1., 1., 100.)
                .unwrap();
            field.add_position_constraint(1, 0., 4.).unwrap();
            field.set_fixed_atoms(&[0]).unwrap();
            let (energy, gradient) = field.energy_gradient().unwrap();
            assert_eq!(energy.to_bits(), expected_energy, "{variant}");
            assert_eq!(gradient[1][0].to_bits(), expected_gradient, "{variant}");
            assert_eq!(gradient[0], [0.; 3]);
            let (converged, _, energy) = field.minimize(20, 1e-4, 1e-6).unwrap();
            assert!(converged);
            assert_eq!(energy.to_bits(), final_energy, "{variant}");
            assert_eq!(field.position(0).unwrap(), [0.; 3]);
            assert_eq!(coordinates.conformers_3d[0].coordinates()[1], [2.5, 0., 0.]);
        }
    }
    #[test]
    fn persistent_angle_constraint_active_capture_and_source_error_order() {
        let initial = vec![[1., 0., 0.], [0.; 3], [0., 1., 0.]];
        let mut field = constraint_only(initial.clone());
        field
            .add_angle_constraint(0, 1, 2, false, 100., 110., 2.)
            .unwrap();
        // Source's 90-degree angle: k*(90-100)^2 = 200.
        assert_eq!(field.energy().unwrap().to_bits(), 200.0f64.to_bits());
        let g = field.gradient().unwrap();
        let derivative = 2.0 * (180.0 / std::f64::consts::PI) * 2.0 * -10.0;
        assert_eq!(
            g,
            [
                [0., -derivative, 0.],
                [derivative, derivative, 0.],
                [-derivative, 0., 0.]
            ]
        );
        let mut relative = constraint_only(initial);
        relative
            .add_angle_constraint(0, 1, 2, true, -10., 0., 2.)
            .unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 0.0f64.to_bits());
        relative.set_position(2, [-1., 0., 0.]).unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 16200.0f64.to_bits());
        let before = relative.energy().unwrap();
        assert_eq!(
            relative
                .add_angle_constraint(9, 1, 2, false, 190., 180., 2.)
                .unwrap_err()
                .kind(),
            MolecularForceFieldErrorKind::InvalidAtomIndex
        );
        assert_eq!(
            relative
                .add_angle_constraint(0, 1, 2, false, 190., 180., 2.)
                .unwrap_err()
                .kind(),
            MolecularForceFieldErrorKind::InvalidConstraintBounds
        );
        let error = relative
            .add_angle_constraint(0, 1, 2, true, 0., 1., 2.)
            .unwrap_err();
        assert_eq!(
            (error.kind(), error.component()),
            (
                MolecularForceFieldErrorKind::InvalidConstraintRange,
                Some(1)
            )
        );
        assert_eq!(relative.energy().unwrap().to_bits(), before.to_bits());
    }
    #[test]
    fn persistent_torsion_constraint_bitexact_wrapped_capture_and_error_order() {
        let rows = vec![[1., 0., 0.], [0.; 3], [0., 1., 0.], [0., 1., 1.]];
        let mut field = constraint_only(rows.clone());
        field
            .add_torsion_constraint(0, 1, 2, 3, false, 0., 0., 2.)
            .unwrap();
        // RDKit 2026.09.1, constraint-only UFF C.C.C.C, captured in ignored
        // source receipt; exact signed dihedral is -90 degrees.
        assert_eq!(field.energy().unwrap().to_bits(), 0x40cfa40000000000);
        assert_eq!(
            field
                .gradient()
                .unwrap()
                .iter()
                .flatten()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            [
                0,
                0,
                0xc0d4249ec28e2462,
                0,
                0,
                0x40d4249ec28e2462,
                0x40d4249ec28e2462,
                0,
                0,
                0xc0d4249ec28e2462,
                0,
                0
            ]
        );
        let mut wrapped = constraint_only(rows.clone());
        wrapped
            .add_torsion_constraint(0, 1, 2, 3, false, 170., 190., 2.)
            .unwrap();
        assert_eq!(wrapped.energy().unwrap().to_bits(), 0x40c9000000000000);
        assert_eq!(
            wrapped.gradient().unwrap()[0][2].to_bits(),
            0x40d1e7a9907e593b
        );
        let mut relative = constraint_only(rows);
        relative
            .add_torsion_constraint(0, 1, 2, 3, true, 0., 0., 2.)
            .unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 0.0f64.to_bits());
        relative.set_position(3, [1., 1., 0.]).unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 16200.0f64.to_bits());
        let before = relative.energy().unwrap();
        assert_eq!(
            relative
                .add_torsion_constraint(8, 1, 2, 3, false, 4., 2., 2.)
                .unwrap_err()
                .kind(),
            MolecularForceFieldErrorKind::InvalidConstraintBounds
        );
        let error = relative
            .add_torsion_constraint(0, 1, 2, 8, false, 0., 0., 2.)
            .unwrap_err();
        assert_eq!(
            (error.kind(), error.atom_index(), error.component()),
            (
                MolecularForceFieldErrorKind::InvalidAtomIndex,
                Some(8),
                Some(3)
            )
        );
        assert_eq!(relative.energy().unwrap().to_bits(), before.to_bits());
    }
    #[test]
    fn persistent_distance_constraint_active_bounds_and_addition_capture() {
        // DistanceConstraint.cpp: res=0.5*k*term*term, endpoint gradient +/-k*term.
        let mut field = constraint_only(vec![[0.; 3], [3., 0., 0.]]);
        field
            .add_distance_constraint(0, 1, false, 1., 2., 4.)
            .unwrap();
        assert_eq!(field.energy().unwrap().to_bits(), 2.0f64.to_bits());
        assert_eq!(field.gradient().unwrap(), [[-4., 0., 0.], [4., 0., 0.]]);
        for distance in [1.0, 1.5, 2.0] {
            field.set_position(1, [distance, 0., 0.]).unwrap();
            assert_eq!(field.energy().unwrap().to_bits(), 0.0f64.to_bits());
        }
        field.set_position(1, [0., 0., 0.]).unwrap();
        assert_eq!(field.energy().unwrap().to_bits(), 2.0f64.to_bits());
        assert_eq!(field.gradient().unwrap(), [[0.; 3]; 2]);
        let mut relative = constraint_only(vec![[0.; 3], [5., 0., 0.]]);
        relative
            .add_distance_constraint(0, 1, true, -2., -1., 4.)
            .unwrap();
        relative.set_position(1, [2., 0., 0.]).unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 2.0f64.to_bits());
        assert_eq!(relative.gradient().unwrap(), [[4., 0., 0.], [-4., 0., 0.]]);
        relative
            .add_distance_constraint(0, 1, true, -5., -4., 2.)
            .unwrap();
        assert_eq!(relative.energy().unwrap().to_bits(), 6.0f64.to_bits());
    }
    #[test]
    fn persistent_distance_constraint_failure_does_not_append_or_rebind() {
        let mut field = constraint_only(vec![[0.; 3], [3., 0., 0.]]);
        field
            .add_distance_constraint(0, 1, false, 1., 2., -4.)
            .unwrap();
        let before = field.energy().unwrap();
        let pointer = field.positions.coordinates().as_ptr();
        let error = field
            .add_distance_constraint(5, 1, false, 4., 2., 3.)
            .unwrap_err();
        assert_eq!(error.kind(), MolecularForceFieldErrorKind::InvalidAtomIndex);
        assert_eq!(
            (error.atom_index(), error.component(), error.expected()),
            (Some(5), Some(0), Some(2))
        );
        assert!(error.source().is_some());
        for bounds in [(4., 2.), (f64::NAN, 2.), (1., f64::NAN)] {
            assert_eq!(
                field
                    .add_distance_constraint(0, 1, false, bounds.0, bounds.1, 3.)
                    .unwrap_err()
                    .kind(),
                MolecularForceFieldErrorKind::InvalidConstraintBounds
            );
        }
        assert!(field.field.as_ref().unwrap().positions().is_empty());
        assert_eq!(field.positions.coordinates().as_ptr(), pointer);
        assert_eq!(field.energy().unwrap().to_bits(), before.to_bits());
        field
            .add_distance_constraint(0, 1, false, 1., 2., 4.)
            .unwrap();
        assert_eq!(field.energy().unwrap().to_bits(), 0.0f64.to_bits());
    }
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
        assert_eq!(converged, expected.outcome.converged());
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
