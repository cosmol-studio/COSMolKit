//! Owned native evaluators; no mirrored or reconstructed force-field state.
use crate::{
    Molecule,
    builder::Points3,
    errors::{error_projection, exception},
    fields::{Field, IndexArray},
    forcefield_params::*,
    group_values::Point3,
    input,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use std::error::Error;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, prelude::*};
crate::enums::enum_value!(MolecularForceFieldErrorKind,"MolecularForceFieldErrorKind"{MissingConformer="missing_conformer",Parameterization="parameterization",InvalidParameterization="invalid_parameterization",InvalidAtomIndex="invalid_atom_index",InvalidFixedAtom="invalid_fixed_atom",CoordinateCount="coordinate_count",CoordinateShape="coordinate_shape",NonFiniteCoordinate="non_finite_coordinate",InvalidTolerance="invalid_tolerance",InvalidConstraintBounds="invalid_constraint_bounds",InvalidConstraintRange="invalid_constraint_range",Preparation="preparation",Construction="construction",Initialization="initialization",Rings="rings",Kernel="kernel"});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "MmffMolPropertiesError | MmffBuilderError | UffParameterError | UffBuilderError | ForceFieldKernelError | ForceFieldError | RingFindingError | OperationError"
    )]
    #[derive(Clone)]
    pub type FieldCause;
}
fn cause<T: Error + 'static>(owner: &ReadonlyValue<T>) -> Option<FieldCause> {
    owner.source().map(|_| {
        let node = owner.child(|e| e.source().expect("immutable force-field source"));
        if node.is::<ck::MmffMolPropertiesError>() {
            return crate::mmff_errors::MmffMolPropertiesError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::MmffMolPropertiesError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        if node.is::<ck::MmffBuilderError>() {
            return crate::forcefield_execution_errors::MmffBuilderError::from_owner(node.child(
                |e| {
                    e.downcast_ref::<ck::MmffBuilderError>()
                        .expect("selected concrete source")
                },
            ))
            .unchecked_into();
        }
        if node.is::<ck::UffParameterError>() {
            return crate::uff_errors::UffParameterError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::UffParameterError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        if node.is::<ck::UffBuilderError>() {
            return crate::uff_errors::UffBuilderError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::UffBuilderError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        if node.is::<ck::ForceFieldKernelError>() {
            return crate::forcefield_kernel_errors::ForceFieldKernelError::from_owner(node.child(
                |e| {
                    e.downcast_ref::<ck::ForceFieldKernelError>()
                        .expect("selected concrete source")
                },
            ))
            .unchecked_into();
        }
        if node.is::<ck::RingFindingError>() {
            return crate::ring_errors::RingFindingError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::RingFindingError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        if node.is::<ck::OperationError>() {
            return crate::operation_errors::OperationError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::OperationError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        if node.is::<ck::ForceFieldError>() {
            return ForceFieldError::from_owner(node.child(|e| {
                e.downcast_ref::<ck::ForceFieldError>()
                    .expect("selected concrete source")
            }))
            .unchecked_into();
        }
        unreachable!("unmodeled native persistent force-field source")
    })
}
error_projection!(
    ForceFieldError,
    ck::ForceFieldError,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> MolecularForceFieldErrorKind {
            self.inner.kind().into()
        }
        #[wasm_bindgen(getter,js_name=requested)]
        pub fn requested(&self) -> Option<usize> {
            self.inner.requested()
        }
        #[wasm_bindgen(getter,js_name=atomIndex)]
        pub fn atom_index(&self) -> Option<usize> {
            self.inner.atom_index()
        }
        #[wasm_bindgen(getter,js_name=component)]
        pub fn component(&self) -> Option<usize> {
            self.inner.component()
        }
        #[wasm_bindgen(getter,js_name=actual)]
        pub fn actual(&self) -> Option<usize> {
            self.inner.actual()
        }
        #[wasm_bindgen(getter,js_name=expected)]
        pub fn expected(&self) -> Option<usize> {
            self.inner.expected()
        }
    },
    source = FieldCause,
    cause = cause::<ck::ForceFieldError>
);
error_projection!(
    MmffForceFieldError,
    ck::MmffForceFieldError,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> MolecularForceFieldErrorKind {
            self.inner.kind().into()
        }
        #[wasm_bindgen(getter,js_name=requested)]
        pub fn requested(&self) -> Option<usize> {
            self.inner.requested()
        }
        #[wasm_bindgen(getter,js_name=atomIndex)]
        pub fn atom_index(&self) -> Option<usize> {
            self.inner.atom_index()
        }
        #[wasm_bindgen(getter,js_name=component)]
        pub fn component(&self) -> Option<usize> {
            self.inner.component()
        }
        #[wasm_bindgen(getter,js_name=actual)]
        pub fn actual(&self) -> Option<usize> {
            self.inner.actual()
        }
        #[wasm_bindgen(getter,js_name=expected)]
        pub fn expected(&self) -> Option<usize> {
            self.inner.expected()
        }
    },
    source = FieldCause,
    cause = cause::<ck::MmffForceFieldError>
);
error_projection!(
    UffForceFieldError,
    ck::UffForceFieldError,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> MolecularForceFieldErrorKind {
            self.inner.kind().into()
        }
        #[wasm_bindgen(getter,js_name=requested)]
        pub fn requested(&self) -> Option<usize> {
            self.inner.requested()
        }
    },
    source = FieldCause,
    cause = cause::<ck::UffForceFieldError>
);

#[wasm_bindgen]
pub struct ForceFieldEnergyGradient(ck::ForceFieldEnergyGradient);
#[wasm_bindgen]
impl ForceFieldEnergyGradient {
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> f64 {
        self.0.energy()
    }
    #[wasm_bindgen(getter)]
    pub fn gradient(&self) -> Result<Ts<Points3>, JsValue> {
        points(self.0.gradient().to_vec())
    }
}
#[wasm_bindgen]
pub struct ForceFieldMinimizeOutcome(ck::ForceFieldMinimizeOutcome);
#[wasm_bindgen]
impl ForceFieldMinimizeOutcome {
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> f64 {
        self.0.energy()
    }
    #[wasm_bindgen(getter)]
    pub fn converged(&self) -> bool {
        self.0.converged()
    }
    #[wasm_bindgen(getter)]
    pub fn iterations(&self) -> u32 {
        self.0.iterations()
    }
}
fn points(value: Vec<[f64; 3]>) -> Result<Ts<Points3>, JsValue> {
    Points3(value.into_iter().map(Point3).collect())
        .into_ts()
        .map_err(|e| wasm_bindgen::JsError::from(e).into())
}
#[wasm_bindgen]
pub struct MolecularForceField(ck::MolecularForceField);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=mmffForceFieldWithParams)]
    pub fn mmff_force_field_with_params(
        &self,
        params: MmffForceFieldParamsInput,
    ) -> Result<MolecularForceField, JsValue> {
        self.mmff_force_field(Some(params))
    }
    #[wasm_bindgen(js_name=uffForceFieldWithParams)]
    pub fn uff_force_field_with_params(
        &self,
        params: UffForceFieldParamsInput,
    ) -> Result<MolecularForceField, JsValue> {
        self.uff_force_field(Some(params))
    }
    #[wasm_bindgen(js_name=mmffForceField)]
    pub fn mmff_force_field(
        &self,
        params: Option<MmffForceFieldParamsInput>,
    ) -> Result<MolecularForceField, JsValue> {
        MmffForceFieldParams::with_input(params.as_ref(), |p| {
            self.0
                .mmff_force_field_with_params(p)
                .map(MolecularForceField)
                .map_err(MmffForceFieldError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=uffForceField)]
    pub fn uff_force_field(
        &self,
        params: Option<UffForceFieldParamsInput>,
    ) -> Result<MolecularForceField, JsValue> {
        UffForceFieldParams::with_input(params.as_ref(), |p| {
            self.0
                .uff_force_field_with_params(p)
                .map(MolecularForceField)
                .map_err(UffForceFieldError::from_rust)
        })
    }
}
#[wasm_bindgen]
impl MolecularForceField {
    #[wasm_bindgen(js_name=minimizeWithParams_)]
    pub fn minimize_with_params_(
        &mut self,
        params: ForceFieldMinimizeParamsInput,
    ) -> Result<ForceFieldMinimizeOutcome, JsValue> {
        self.minimize_(Some(params))
    }
    #[wasm_bindgen]
    pub fn position(&self, atom: Ts<input::IndexInput>) -> Result<Ts<Point3>, JsValue> {
        let v = self
            .0
            .position(ck::AtomId::new(input::read::<usize>(atom)?))
            .map_err(ForceFieldError::from_rust)?;
        Point3(v)
            .into_ts()
            .map_err(|e| wasm_bindgen::JsError::from(e).into())
    }
    #[wasm_bindgen(js_name=setPosition_)]
    pub fn set_position_(
        &mut self,
        atom: Ts<input::IndexInput>,
        position: Ts<Point3>,
    ) -> Result<(), JsValue> {
        self.0
            .set_position_(
                ck::AtomId::new(input::read::<usize>(atom)?),
                input::decode(position)?.0,
            )
            .map_err(ForceFieldError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn positions(&self) -> Result<Ts<Points3>, JsValue> {
        points(self.0.positions())
    }
    #[wasm_bindgen(js_name=setPositions_)]
    pub fn set_positions_(&mut self, positions: Ts<Points3>) -> Result<(), JsValue> {
        let p = input::decode(positions)?
            .0
            .into_iter()
            .map(|p| p.0)
            .collect::<Vec<_>>();
        self.0
            .set_positions_(&p)
            .map_err(ForceFieldError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=fixedAtoms)]
    pub fn fixed_atoms(&self) -> Vec<usize> {
        self.0
            .fixed_atoms()
            .into_iter()
            .map(|id| id.index())
            .collect()
    }
    #[wasm_bindgen(js_name=setFixedAtoms_)]
    pub fn set_fixed_atoms_(&mut self, atoms: Ts<IndexArray>) -> Result<(), JsValue> {
        let atoms = Vec::<usize>::input(atoms)?
            .into_iter()
            .map(ck::AtomId::new)
            .collect::<Vec<_>>();
        self.0
            .set_fixed_atoms_(&atoms)
            .map_err(ForceFieldError::from_rust)
    }
    #[wasm_bindgen]
    pub fn energy(&self) -> Result<f64, JsValue> {
        self.0.energy().map_err(ForceFieldError::from_rust)
    }
    #[wasm_bindgen]
    pub fn gradient(&self) -> Result<Ts<Points3>, JsValue> {
        points(self.0.gradient().map_err(ForceFieldError::from_rust)?)
    }
    #[wasm_bindgen(js_name=gradientUnconstrained)]
    pub fn gradient_unconstrained(&self) -> Result<Ts<Points3>, JsValue> {
        points(
            self.0
                .gradient_unconstrained()
                .map_err(ForceFieldError::from_rust)?,
        )
    }
    #[wasm_bindgen(js_name=energyGradient)]
    pub fn energy_gradient(&self) -> Result<ForceFieldEnergyGradient, JsValue> {
        self.0
            .energy_gradient()
            .map(ForceFieldEnergyGradient)
            .map_err(ForceFieldError::from_rust)
    }
    #[wasm_bindgen(js_name=minimize_)]
    pub fn minimize_(
        &mut self,
        params: Option<ForceFieldMinimizeParamsInput>,
    ) -> Result<ForceFieldMinimizeOutcome, JsValue> {
        ForceFieldMinimizeParams::with_input(params.as_ref(), |p| {
            self.0
                .minimize_with_params_(p)
                .map(ForceFieldMinimizeOutcome)
                .map_err(ForceFieldError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addDistanceConstraint_)]
    pub fn add_distance_constraint_(
        &mut self,
        atom1: Ts<input::IndexInput>,
        atom2: Ts<input::IndexInput>,
        params: ForceFieldDistanceConstraintParamsInput,
    ) -> Result<(), JsValue> {
        let atom1 = ck::AtomId::new(input::read::<usize>(atom1)?);
        let atom2 = ck::AtomId::new(input::read::<usize>(atom2)?);
        ForceFieldDistanceConstraintParams::with_input(Some(&params), |p| {
            self.0
                .add_distance_constraint_(atom1, atom2, p)
                .map_err(ForceFieldError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addAngleConstraint_)]
    pub fn add_angle_constraint_(
        &mut self,
        atom1: Ts<input::IndexInput>,
        atom2: Ts<input::IndexInput>,
        atom3: Ts<input::IndexInput>,
        params: ForceFieldAngleConstraintParamsInput,
    ) -> Result<(), JsValue> {
        let atom1 = ck::AtomId::new(input::read::<usize>(atom1)?);
        let atom2 = ck::AtomId::new(input::read::<usize>(atom2)?);
        let atom3 = ck::AtomId::new(input::read::<usize>(atom3)?);
        ForceFieldAngleConstraintParams::with_input(Some(&params), |p| {
            self.0
                .add_angle_constraint_(atom1, atom2, atom3, p)
                .map_err(ForceFieldError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addTorsionConstraint_)]
    pub fn add_torsion_constraint_(
        &mut self,
        atom1: Ts<input::IndexInput>,
        atom2: Ts<input::IndexInput>,
        atom3: Ts<input::IndexInput>,
        atom4: Ts<input::IndexInput>,
        params: ForceFieldTorsionConstraintParamsInput,
    ) -> Result<(), JsValue> {
        let atom1 = ck::AtomId::new(input::read::<usize>(atom1)?);
        let atom2 = ck::AtomId::new(input::read::<usize>(atom2)?);
        let atom3 = ck::AtomId::new(input::read::<usize>(atom3)?);
        let atom4 = ck::AtomId::new(input::read::<usize>(atom4)?);
        ForceFieldTorsionConstraintParams::with_input(Some(&params), |p| {
            self.0
                .add_torsion_constraint_(atom1, atom2, atom3, atom4, p)
                .map_err(ForceFieldError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addPositionConstraint_)]
    pub fn add_position_constraint_(
        &mut self,
        atom: Ts<input::IndexInput>,
        params: ForceFieldPositionConstraintParamsInput,
    ) -> Result<(), JsValue> {
        let atom = ck::AtomId::new(input::read::<usize>(atom)?);
        ForceFieldPositionConstraintParams::with_input(Some(&params), |p| {
            self.0
                .add_position_constraint_(atom, p)
                .map_err(ForceFieldError::from_rust)
        })
    }
}
