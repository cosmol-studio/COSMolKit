//! Native force-field properties, diagnostics and value-returning optimization.
use crate::{
    Molecule,
    errors::child,
    errors::{enum_error, exception},
    forcefield_params::*,
    input,
    mmff_errors::MmffMolPropertiesError,
    operation_errors::OperationError,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
ck::__ck_configuration_UffOptimizationParams!(
    crate::project_parameters,
    ck::UffOptimizationParams,
    {},
    { UffOptimizationParams }
);
ck::__ck_configuration_UffConformerOptimizationParams!(
    crate::project_parameters,
    ck::UffConformerOptimizationParams,
    {},
    { UffConformerOptimizationParams }
);
ck::__ck_configuration_UffEvaluationParams!(
    crate::project_parameters,
    ck::UffEvaluationParams,
    {},
    { UffEvaluationParams }
);
ck::__ck_configuration_MmffOptimizationParams!(
    crate::project_parameters,
    ck::MmffOptimizationParams,
    {},
    { MmffOptimizationParams }
);
ck::__ck_configuration_MmffConformerOptimizationParams!(
    crate::project_parameters,
    ck::MmffConformerOptimizationParams,
    {},
    { MmffConformerOptimizationParams }
);
ck::__ck_configuration_MmffEvaluationParams!(
    crate::project_parameters,
    ck::MmffEvaluationParams,
    {},
    { MmffEvaluationParams }
);

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "OperationError | UffParameterError")]
    #[derive(Clone)]
    pub type QueryCause;
}
enum_error!(UffParameterQueryError,UffParameterQueryErrorKind,ck::UffParameterQueryError,{Cache=>ck::UffParameterQueryError::Cache(_),Parameters=>ck::UffParameterQueryError::Parameters(_)},{},source=QueryCause,cause=|o:&ReadonlyValue<ck::UffParameterQueryError>|Some(match o.as_ref(){ck::UffParameterQueryError::Cache(_)=>child!(o,ck::UffParameterQueryError::Cache,OperationError),ck::UffParameterQueryError::Parameters(_)=>child!(o,ck::UffParameterQueryError::Parameters,crate::uff_errors::UffParameterError)}));
#[wasm_bindgen]
pub struct MmffAtomProperties(ck::MmffAtomProperties);
#[wasm_bindgen]
impl MmffAtomProperties {
    #[wasm_bindgen(constructor)]
    pub fn new(
        atom_type: Ts<input::U8Input>,
        formal_charge: Ts<input::F64Input>,
        partial_charge: Ts<input::F64Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MmffAtomProperties {
            atom_type: input::read::<u8>(atom_type)?,
            formal_charge: input::read::<f64>(formal_charge)?,
            partial_charge: input::read::<f64>(partial_charge)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=atomType)]
    pub fn atom_type(&self) -> u8 {
        self.0.atom_type
    }
    #[wasm_bindgen(getter,js_name=formalCharge)]
    pub fn formal_charge(&self) -> f64 {
        self.0.formal_charge
    }
    #[wasm_bindgen(getter,js_name=partialCharge)]
    pub fn partial_charge(&self) -> f64 {
        self.0.partial_charge
    }
}
#[wasm_bindgen]
pub struct MmffProperties(ck::MmffProperties);
#[wasm_bindgen]
impl MmffProperties {
    #[wasm_bindgen(getter,js_name=isValid)]
    pub fn is_valid(&self) -> bool {
        self.0.is_valid()
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> MmffVariant {
        self.0.variant().into()
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Vec<MmffAtomProperties> {
        self.0
            .atoms()
            .iter()
            .cloned()
            .map(MmffAtomProperties)
            .collect()
    }
    #[wasm_bindgen(js_name=atomType)]
    pub fn atom_type(&self, index: Ts<input::IndexInput>) -> Result<u8, JsValue> {
        self.0
            .atom_type(input::read::<usize>(index)?)
            .map_err(MmffMolPropertiesError::from_rust)
    }
    #[wasm_bindgen(js_name=formalCharge)]
    pub fn formal_charge(&self, index: Ts<input::IndexInput>) -> Result<f64, JsValue> {
        self.0
            .formal_charge(input::read::<usize>(index)?)
            .map_err(MmffMolPropertiesError::from_rust)
    }
    #[wasm_bindgen(js_name=partialCharge)]
    pub fn partial_charge(&self, index: Ts<input::IndexInput>) -> Result<f64, JsValue> {
        self.0
            .partial_charge(input::read::<usize>(index)?)
            .map_err(MmffMolPropertiesError::from_rust)
    }
}
#[wasm_bindgen]
pub struct MmffEnergyGradient(ck::MmffEnergyGradient);
#[wasm_bindgen]
impl MmffEnergyGradient {
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> f64 {
        self.0.energy()
    }
    #[wasm_bindgen(getter)]
    pub fn gradient(&self) -> Vec<f64> {
        self.0.gradient().to_vec()
    }
}
#[wasm_bindgen]
pub struct UffEnergyGradient(ck::UffEnergyGradient);
#[wasm_bindgen]
impl UffEnergyGradient {
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> f64 {
        self.0.energy()
    }
    #[wasm_bindgen(getter)]
    pub fn gradient(&self) -> Vec<f64> {
        self.0.gradient().to_vec()
    }
}
crate::enums::enum_value!(OptimizationStatus,"OptimizationStatus"{Converged="converged",NotConverged="not_converged",MissingParameters="missing_parameters"});
#[wasm_bindgen]
pub struct OptimizationOutcome(ck::OptimizationOutcome);
#[wasm_bindgen]
impl OptimizationOutcome {
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> OptimizationStatus {
        self.0.status().into()
    }
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> Option<f64> {
        self.0.energy()
    }
    #[wasm_bindgen(getter)]
    pub fn converged(&self) -> bool {
        self.0.converged()
    }
    #[wasm_bindgen(getter,js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.0.needs_more()
    }
    #[wasm_bindgen(getter,js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.0.status_code()
    }
}
#[wasm_bindgen]
pub struct UffOptimizationResult(ck::UffOptimizationResult);
#[wasm_bindgen]
impl UffOptimizationResult {
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome(self.0.outcome)
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeResult(ck::MmffOptimizeMoleculeResult);
#[wasm_bindgen]
impl MmffOptimizeMoleculeResult {
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome(self.0.outcome)
    }
}
#[wasm_bindgen]
pub struct UffConformerResult(ck::UffConformerResult);
#[wasm_bindgen]
impl UffConformerResult {
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> usize {
        self.0.conformer_id
    }
    #[wasm_bindgen(getter)]
    pub fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome(self.0.outcome)
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeConfResult(ck::MmffOptimizeMoleculeConfResult);
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfResult {
    #[wasm_bindgen(getter)]
    pub fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome(self.0.outcome)
    }
}
#[wasm_bindgen]
pub struct UffConformerOptimizationResult(ck::UffConformerOptimizationResult);
#[wasm_bindgen]
impl UffConformerOptimizationResult {
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule.clone())
    }
    #[wasm_bindgen(getter,js_name=conformers)]
    pub fn conformers(&self) -> Vec<UffConformerResult> {
        self.0
            .conformers
            .iter()
            .copied()
            .map(UffConformerResult)
            .collect()
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeConfsResult(ck::MmffOptimizeMoleculeConfsResult);
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfsResult {
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule.clone())
    }
    #[wasm_bindgen(getter,js_name=conformerResults)]
    pub fn conformer_results(&self) -> Vec<MmffOptimizeMoleculeConfResult> {
        self.0
            .conformer_results
            .iter()
            .copied()
            .map(MmffOptimizeMoleculeConfResult)
            .collect()
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=uffHasAllMoleculeParams)]
    pub fn uff_has_all_molecule_params(&self) -> Result<bool, JsValue> {
        self.0
            .uff_has_all_molecule_params()
            .map_err(UffParameterQueryError::from_rust)
    }
    #[wasm_bindgen(js_name=mmffHasAllMoleculeParams)]
    pub fn mmff_has_all_molecule_params(&self) -> Result<bool, JsValue> {
        self.0
            .mmff_has_all_molecule_params()
            .map_err(MmffMolPropertiesError::from_rust)
    }
    #[wasm_bindgen(js_name=mmffProperties)]
    pub fn mmff_properties(
        &self,
        params: Option<MmffPropertiesParamsInput>,
    ) -> Result<MmffProperties, JsValue> {
        MmffPropertiesParams::with_input(params.as_ref(), |p| {
            self.0
                .mmff_properties_with_params(p)
                .map(MmffProperties)
                .map_err(MmffMolPropertiesError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=mmffEnergyGradient)]
    pub fn mmff_energy_gradient(
        &self,
        params: Option<MmffEvaluationParamsInput>,
    ) -> Result<Option<MmffEnergyGradient>, JsValue> {
        MmffEvaluationParams::with_input(params.as_ref(), |p| {
            self.0
                .mmff_energy_gradient_with_params(p)
                .map(|v| v.map(MmffEnergyGradient))
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withUffOptimized)]
    pub fn with_uff_optimized(
        &self,
        params: Option<UffOptimizationParamsInput>,
    ) -> Result<UffOptimizationResult, JsValue> {
        UffOptimizationParams::with_input(params.as_ref(), |p| {
            self.0
                .with_uff_optimized_with_params(p)
                .map(UffOptimizationResult)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withUffOptimizedConformers)]
    pub fn with_uff_optimized_conformers(
        &self,
        params: Option<UffConformerOptimizationParamsInput>,
    ) -> Result<UffConformerOptimizationResult, JsValue> {
        UffConformerOptimizationParams::with_input(params.as_ref(), |p| {
            self.0
                .with_uff_optimized_conformers_with_params(p)
                .map(UffConformerOptimizationResult)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withMmffOptimized)]
    pub fn with_mmff_optimized(
        &self,
        params: Option<MmffOptimizationParamsInput>,
    ) -> Result<MmffOptimizeMoleculeResult, JsValue> {
        MmffOptimizationParams::with_input(params.as_ref(), |p| {
            self.0
                .with_mmff_optimized_with_params(p)
                .map(MmffOptimizeMoleculeResult)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withMmffOptimizedConformers)]
    pub fn with_mmff_optimized_conformers(
        &self,
        params: Option<MmffConformerOptimizationParamsInput>,
    ) -> Result<MmffOptimizeMoleculeConfsResult, JsValue> {
        MmffConformerOptimizationParams::with_input(params.as_ref(), |p| {
            self.0
                .with_mmff_optimized_conformers_with_params(p)
                .map(MmffOptimizeMoleculeConfsResult)
                .map_err(OperationError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=uffEnergyGradient)]
    pub fn uff_energy_gradient(
        &self,
        params: Option<UffEvaluationParamsInput>,
    ) -> Result<UffEnergyGradient, JsValue> {
        UffEvaluationParams::with_input(params.as_ref(), |p| {
            self.0
                .uff_energy_gradient_with_params(p)
                .map(UffEnergyGradient)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=uffEnergyGradientWithParams)]
    pub fn uff_energy_gradient_with_params(
        &self,
        params: UffEvaluationParamsInput,
    ) -> Result<UffEnergyGradient, JsValue> {
        self.uff_energy_gradient(Some(params))
    }
    #[wasm_bindgen(js_name=mmffEnergyGradientWithParams)]
    pub fn mmff_energy_gradient_with_params(
        &self,
        params: MmffEvaluationParamsInput,
    ) -> Result<Option<MmffEnergyGradient>, JsValue> {
        self.mmff_energy_gradient(Some(params))
    }
    #[wasm_bindgen(js_name=mmffPropertiesWithParams)]
    pub fn mmff_properties_with_params(
        &self,
        params: MmffPropertiesParamsInput,
    ) -> Result<MmffProperties, JsValue> {
        self.mmff_properties(Some(params))
    }
    #[wasm_bindgen(js_name=withUffOptimizedWithParams)]
    pub fn with_uff_optimized_with_params(
        &self,
        params: UffOptimizationParamsInput,
    ) -> Result<UffOptimizationResult, JsValue> {
        self.with_uff_optimized(Some(params))
    }
    #[wasm_bindgen(js_name=withUffOptimizedConformersWithParams)]
    pub fn with_uff_optimized_conformers_with_params(
        &self,
        params: UffConformerOptimizationParamsInput,
    ) -> Result<UffConformerOptimizationResult, JsValue> {
        self.with_uff_optimized_conformers(Some(params))
    }
    #[wasm_bindgen(js_name=withMmffOptimizedWithParams)]
    pub fn with_mmff_optimized_with_params(
        &self,
        params: MmffOptimizationParamsInput,
    ) -> Result<MmffOptimizeMoleculeResult, JsValue> {
        self.with_mmff_optimized(Some(params))
    }
    #[wasm_bindgen(js_name=withMmffOptimizedConformersWithParams)]
    pub fn with_mmff_optimized_conformers_with_params(
        &self,
        params: MmffConformerOptimizationParamsInput,
    ) -> Result<MmffOptimizeMoleculeConfsResult, JsValue> {
        self.with_mmff_optimized_conformers(Some(params))
    }
}

// Detached result constructors preserve the same native values as their producers.
use crate::fields::Field;
crate::fields::value_field!(
    ck::OptimizationOutcome,
    OptimizationOutcome,
    OptimizationOutcomeInput,
    "OptimizationOutcome"
);
crate::fields::value_field!(
    ck::UffConformerResult,
    UffConformerResult,
    UffConformerResultInput,
    "UffConformerResult"
);
crate::fields::value_field!(
    ck::MmffOptimizeMoleculeConfResult,
    MmffOptimizeMoleculeConfResult,
    MmffOptimizeMoleculeConfResultInput,
    "MmffOptimizeMoleculeConfResult"
);
crate::fields::value_array!(
    ck::UffConformerResult,
    UffConformerResultArray,
    "ReadonlyArray<UffConformerResult>"
);
crate::fields::value_array!(
    ck::MmffOptimizeMoleculeConfResult,
    MmffOptimizeMoleculeConfResultArray,
    "ReadonlyArray<MmffOptimizeMoleculeConfResult>"
);
#[wasm_bindgen]
impl OptimizationOutcome {
    #[wasm_bindgen(js_name=convergedWithEnergy)]
    pub fn converged_with_energy(energy: Ts<input::F64Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::OptimizationOutcome::converged_with_energy(
            input::read::<f64>(energy)?,
        )))
    }
    #[wasm_bindgen(js_name=notConvergedWithEnergy)]
    pub fn not_converged_with_energy(energy: Ts<input::F64Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::OptimizationOutcome::not_converged_with_energy(
            input::read::<f64>(energy)?,
        )))
    }
    #[wasm_bindgen(js_name=missingParameters)]
    pub fn missing_parameters() -> Self {
        Self(ck::OptimizationOutcome::missing_parameters())
    }
}
#[wasm_bindgen]
impl UffOptimizationResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        molecule: crate::MoleculeInput,
        outcome: OptimizationOutcomeInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::UffOptimizationResult {
            molecule: ck::Molecule::input(molecule)?,
            outcome: ck::OptimizationOutcome::input(outcome)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.0.status_code()
    }
    #[wasm_bindgen(getter,js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.0.needs_more()
    }
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> Option<f64> {
        self.0.energy()
    }
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        molecule: crate::MoleculeInput,
        outcome: OptimizationOutcomeInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MmffOptimizeMoleculeResult {
            molecule: ck::Molecule::input(molecule)?,
            outcome: ck::OptimizationOutcome::input(outcome)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.0.status_code()
    }
    #[wasm_bindgen(getter,js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.0.needs_more()
    }
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> Option<f64> {
        self.0.energy()
    }
}
#[wasm_bindgen]
impl UffConformerResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        conformer_id: Ts<input::IndexInput>,
        outcome: OptimizationOutcomeInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::UffConformerResult {
            conformer_id: input::read::<usize>(conformer_id)?,
            outcome: ck::OptimizationOutcome::input(outcome)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.0.status_code()
    }
    #[wasm_bindgen(getter,js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.0.needs_more()
    }
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> Option<f64> {
        self.0.energy()
    }
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfResult {
    #[wasm_bindgen(constructor)]
    pub fn new(outcome: OptimizationOutcomeInput) -> Result<Self, JsValue> {
        Ok(Self(ck::MmffOptimizeMoleculeConfResult {
            outcome: ck::OptimizationOutcome::input(outcome)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.0.status_code()
    }
    #[wasm_bindgen(getter,js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.0.needs_more()
    }
    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> Option<f64> {
        self.0.energy()
    }
}
#[wasm_bindgen]
impl UffConformerOptimizationResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        molecule: crate::MoleculeInput,
        values: UffConformerResultArray,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::UffConformerOptimizationResult {
            molecule: ck::Molecule::input(molecule)?,
            conformers: Vec::<ck::UffConformerResult>::input(values)?,
        }))
    }
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfsResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        molecule: crate::MoleculeInput,
        values: MmffOptimizeMoleculeConfResultArray,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MmffOptimizeMoleculeConfsResult {
            molecule: ck::Molecule::input(molecule)?,
            conformer_results: Vec::<ck::MmffOptimizeMoleculeConfResult>::input(values)?,
        }))
    }
}
#[wasm_bindgen]
impl UffConformerOptimizationResult {
    #[wasm_bindgen(getter,js_name=conformerResults)]
    pub fn conformer_results(&self) -> Vec<UffConformerResult> {
        self.0
            .conformer_results()
            .iter()
            .copied()
            .map(UffConformerResult)
            .collect()
    }
}
#[wasm_bindgen]
impl UffEnergyGradient {
    #[wasm_bindgen(constructor)]
    pub fn new(
        energy: Ts<input::F64Input>,
        gradient: Ts<crate::fields::FloatArray>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::UffEnergyGradient {
            energy: input::read::<f64>(energy)?,
            gradient: Vec::<f64>::input(gradient)?,
        }))
    }
}
#[wasm_bindgen]
impl MmffEnergyGradient {
    #[wasm_bindgen(constructor)]
    pub fn new(
        energy: Ts<input::F64Input>,
        gradient: Ts<crate::fields::FloatArray>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MmffEnergyGradient {
            energy: input::read::<f64>(energy)?,
            gradient: Vec::<f64>::input(gradient)?,
        }))
    }
}
