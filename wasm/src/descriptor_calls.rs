//! Descriptor queries delegate once and retain native cache/force behavior.
use crate::{Molecule, descriptor_errors::DescriptorReadError, fields::Field, input};
use wasm_bindgen::prelude::*;
#[derive(serde::Serialize, tsify::Tsify)]
#[serde(transparent)]
pub struct HallKierContributions((f64, Vec<f64>));
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=hallKierAlphaWithContributions)]
    pub fn hall_kier_alpha_with_contributions(
        &self,
    ) -> Result<tsify::Ts<HallKierContributions>, JsValue> {
        use tsify::Tsify;
        HallKierContributions(
            self.0
                .hall_kier_alpha_with_contributions()
                .map_err(DescriptorReadError::from_rust)?,
        )
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct CrippenTotals(cosmolkit::CrippenTotals);
#[wasm_bindgen]
impl CrippenTotals {
    #[wasm_bindgen(constructor)]
    pub fn new(
        logp: tsify::Ts<input::F64Input>,
        molar_refractivity: tsify::Ts<input::F64Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(cosmolkit::CrippenTotals {
            logp: input::read::<f64>(logp)?,
            molar_refractivity: input::read::<f64>(molar_refractivity)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn logp(&self) -> f64 {
        self.0.logp
    }
    #[wasm_bindgen(getter,js_name=molarRefractivity)]
    pub fn molar_refractivity(&self) -> f64 {
        self.0.molar_refractivity
    }
}
#[wasm_bindgen]
pub struct LabuteAsaContributions(cosmolkit::LabuteAsaContributions);
#[wasm_bindgen]
impl LabuteAsaContributions {
    #[wasm_bindgen(constructor)]
    pub fn new(
        asa: tsify::Ts<input::F64Input>,
        atom_contributions: tsify::Ts<crate::fields::FloatArray>,
        hydrogen_contribution: tsify::Ts<input::F64Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(cosmolkit::LabuteAsaContributions {
            asa: input::read::<f64>(asa)?,
            atom_contributions: Vec::<f64>::input(atom_contributions)?,
            hydrogen_contribution: input::read::<f64>(hydrogen_contribution)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn asa(&self) -> f64 {
        self.0.asa
    }
    #[wasm_bindgen(getter,js_name=atomContributions)]
    pub fn atom_contributions(&self) -> Vec<f64> {
        self.0.atom_contributions.clone()
    }
    #[wasm_bindgen(getter,js_name=hydrogenContribution)]
    pub fn hydrogen_contribution(&self) -> f64 {
        self.0.hydrogen_contribution
    }
}
// Use the existing explicit adapter hook for native checked accessors.
macro_rules! project_charge_params {
    (@resume {} @adapter $(#[$doc:meta])* $native:path {
        $(#[$count_doc:meta])* $count:ident: $count_ty:ty,
        $(#[$strict_doc:meta])* $strict:ident: $strict_ty:ty
    } methods {}) => {
        paste::paste! {
            #[derive(serde::Deserialize, tsify::Tsify)]
            #[serde(default, rename_all = "camelCase")]
            pub struct GasteigerChargeOptions {
                #[tsify(optional)]
                $count: i32,
                #[tsify(optional, type = "boolean")]
                $strict: $strict_ty,
            }

            impl Default for GasteigerChargeOptions {
                fn default() -> Self {
                    let params = cosmolkit::GasteigerChargeParams::default();
                    Self {
                        // Native checked storage guarantees the signed-int range.
                        $count: params.$count() as i32,
                        $strict: params.$strict(),
                    }
                }
            }

            #[wasm_bindgen]
            extern "C" {
                #[wasm_bindgen(typescript_type = "GasteigerChargeParams | GasteigerChargeOptions")]
                pub type GasteigerConfigInput;

                #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
                fn borrow_charge_params(
                    visit: &mut dyn FnMut(&GasteigerChargeParams),
                    receiver: &JsValue,
                    args: &js_sys::Array,
                ) -> Result<JsValue, JsValue>;
            }


            impl GasteigerChargeOptions {
                fn into_params(self) -> Result<cosmolkit::GasteigerChargeParams, JsValue> {
                    let iterations = i64::from(self.$count);
                    cosmolkit::GasteigerChargeParams::new(iterations, self.$strict)
                        .map_err(|error| js_sys::RangeError::new(&error.to_string()).into())
                }
            }

            impl GasteigerConfigInput {
                fn params(&self) -> Result<cosmolkit::GasteigerChargeParams, JsValue> {
                    let raw = self.as_ref();
                    if input::is_plain_options(raw) {
                        input::check_options(raw, &["iterations", "throwOnParamFailure"])?;
                        let mut fields = GasteigerChargeOptions::default();
                        let iterations = js_sys::Reflect::get(raw, &"iterations".into())?;
                        if !iterations.is_undefined() {
                            fields.$count = i32::decode(&iterations)?;
                        }
                        let strict = js_sys::Reflect::get(raw, &"throwOnParamFailure".into())?;
                        if !strict.is_undefined() {
                            fields.$strict = bool::decode(&strict)?;
                        }
                        return fields.into_params();
                    }
                    input::with_instance(
                        raw,
                        "GasteigerChargeParams",
                        |visit, receiver, args| borrow_charge_params(&mut |p| visit(p), receiver, args),
                        |p| Ok(p.0),
                    )
                }
            }

            #[wasm_bindgen]
            pub struct GasteigerChargeParams(cosmolkit::GasteigerChargeParams);

            #[wasm_bindgen]
            impl GasteigerChargeParams {
                #[wasm_bindgen(constructor)]
                pub fn new(
                     $count: Option<tsify::Ts<input::I32Input>>,
                     $strict: Option<tsify::Ts<input::BoolInput>>,
                ) -> Result<Self, JsValue> {
                    let defaults = GasteigerChargeOptions::default();
                    GasteigerChargeOptions {
                        $count: $count
                            .map(input::read::<i32>)
                            .transpose()?
                            .unwrap_or(defaults.$count),
                        $strict: $strict
                            .map(input::read::<bool>)
                            .transpose()?
                            .unwrap_or(defaults.$strict),
                    }
                    .into_params()
                    .map(Self)
                }

                #[wasm_bindgen(getter)]
                $(#[$count_doc])*
                pub fn $count(&self) -> $count_ty {
                    self.0.$count()
                }

                #[wasm_bindgen(setter)]
                pub fn [<set_ $count>](&mut self, value: tsify::Ts<input::I32Input>) -> Result<(), JsValue> {
                    let value = i64::from(input::read::<i32>(value)?);
                    self.0
                        .[<set_ $count>](value)
                        .map_err(|error| js_sys::RangeError::new(&error.to_string()).into())
                }

                #[wasm_bindgen(getter, js_name = throwOnParamFailure)]
                $(#[$strict_doc])*
                pub fn $strict(&self) -> $strict_ty {
                    self.0.$strict()
                }

                #[wasm_bindgen(setter, js_name = throwOnParamFailure)]
                pub fn [<set_ $strict>](
                    &mut self,
                    value: tsify::Ts<input::BoolInput>,
                ) -> Result<(), JsValue> {
                    self.0
                        .[<set_ $strict>](input::read::<bool>(value)?);
                    Ok(())
                }
            }

        }
    };
}
cosmolkit::__ck_configuration_GasteigerChargeParams!(
    project_charge_params,
    cosmolkit::GasteigerChargeParams,
    {},
    {}
);

#[wasm_bindgen]
pub struct GasteigerCharge(cosmolkit::GasteigerCharge);

#[wasm_bindgen]
impl GasteigerCharge {
    #[wasm_bindgen(getter, js_name = atomCharge)]
    pub fn atom_charge(&self) -> f64 {
        self.0.atom_charge()
    }

    #[wasm_bindgen(getter, js_name = attachedHydrogenCharge)]
    pub fn attached_hydrogen_charge(&self) -> f64 {
        self.0.attached_hydrogen_charge()
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=gasteigerCharges)]
    pub fn gasteiger_charges(
        &self,
        options: Option<GasteigerConfigInput>,
    ) -> Result<Vec<GasteigerCharge>, JsValue> {
        let params = options
            .as_ref()
            .map(GasteigerConfigInput::params)
            .transpose()?
            .unwrap_or_default();
        self.gasteiger_charges_with_params(&GasteigerChargeParams(params))
    }

    #[wasm_bindgen(js_name=gasteigerChargesWithParams)]
    pub fn gasteiger_charges_with_params(
        &self,
        params: &GasteigerChargeParams,
    ) -> Result<Vec<GasteigerCharge>, JsValue> {
        self.0
            .gasteiger_charges_with_params(&params.0)
            .map(|rows| rows.into_iter().map(GasteigerCharge).collect())
            .map_err(DescriptorReadError::from_rust)
    }

    #[wasm_bindgen(js_name=peoeVsa)]
    pub fn peoe_vsa(&self) -> Result<Vec<f64>, JsValue> {
        self.0.peoe_vsa().map_err(DescriptorReadError::from_rust)
    }

    #[wasm_bindgen(js_name=numRotatableBonds)]
    pub fn num_rotatable_bonds(
        &self,
        options: Option<crate::enums::RotatableBondsOptionsInput>,
    ) -> Result<u32, JsValue> {
        let options = options
            .map(cosmolkit::RotatableBondsOptions::input)
            .transpose()?
            .unwrap_or_default();
        self.0
            .num_rotatable_bonds_with_params(&options)
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=crippenDescriptors)]
    pub fn crippen_descriptors(
        &self,
        include_hydrogens: Option<tsify::Ts<input::BoolInput>>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<CrippenTotals, JsValue> {
        self.0
            .crippen_descriptors_with_params(
                include_hydrogens
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(true),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map(CrippenTotals)
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=labuteAsaContributions)]
    pub fn labute_asa_contributions(
        &self,
        include_hydrogens: Option<tsify::Ts<input::BoolInput>>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<LabuteAsaContributions, JsValue> {
        self.0
            .labute_asa_contributions_with_params(
                include_hydrogens
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(true),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map(LabuteAsaContributions)
            .map_err(DescriptorReadError::from_rust)
    }
}
macro_rules! queries{($($name:ident,$js:ident,$out:ty);+$(;)?)=>{#[wasm_bindgen]impl Molecule{$(
    #[wasm_bindgen(js_name=$js)]pub fn $name(&self)->Result<$out,JsValue>{self.0.$name().map_err(DescriptorReadError::from_rust)}
)+}};}
queries! {
 bcut2d,bcut2d,Vec<f64>;
 estate_indices,estateIndices,Vec<f64>;
 max_estate_index,maxEstateIndex,f64;
 min_estate_index,minEstateIndex,f64;
 max_abs_estate_index,maxAbsEstateIndex,f64;
 min_abs_estate_index,minAbsEstateIndex,f64;
 estate_vsa,estateVsa,Vec<f64>;
 vsa_estate,vsaEstate,Vec<f64>;
 estate_vsa_1,estateVsa1,f64;
 estate_vsa_2,estateVsa2,f64;
 estate_vsa_3,estateVsa3,f64;
 estate_vsa_4,estateVsa4,f64;
 estate_vsa_5,estateVsa5,f64;
 estate_vsa_6,estateVsa6,f64;
 estate_vsa_7,estateVsa7,f64;
 estate_vsa_8,estateVsa8,f64;
 estate_vsa_9,estateVsa9,f64;
 estate_vsa_10,estateVsa10,f64;
 estate_vsa_11,estateVsa11,f64;
 vsa_estate_1,vsaEstate1,f64;
 vsa_estate_2,vsaEstate2,f64;
 vsa_estate_3,vsaEstate3,f64;
 vsa_estate_4,vsaEstate4,f64;
 vsa_estate_5,vsaEstate5,f64;
 vsa_estate_6,vsaEstate6,f64;
 vsa_estate_7,vsaEstate7,f64;
 vsa_estate_8,vsaEstate8,f64;
 vsa_estate_9,vsaEstate9,f64;
 vsa_estate_10,vsaEstate10,f64;
 chi_0,chi0,f64;chi_1,chi1,f64;hall_kier_alpha,hallKierAlpha,f64;
 kappa_1,kappa1,f64;kappa_2,kappa2,f64;kappa_3,kappa3,f64;phi,phi,f64;
 num_heavy_atoms,numHeavyAtoms,u32;total_atom_count,totalAtomCount,u32;
 lipinski_hba,lipinskiHba,u32;lipinski_hbd,lipinskiHbd,u32;fraction_csp3,fractionCsp3,f64;
 num_heteroatoms,numHeteroatoms,u32;num_hba,numHba,u32;num_hbd,numHbd,u32;
 num_heterocycles,numHeterocycles,u32;num_aromatic_rings,numAromaticRings,u32;
 num_saturated_rings,numSaturatedRings,u32;num_aliphatic_rings,numAliphaticRings,u32;
 num_aromatic_heterocycles,numAromaticHeterocycles,u32;num_aromatic_carbocycles,numAromaticCarbocycles,u32;
 num_aliphatic_heterocycles,numAliphaticHeterocycles,u32;num_aliphatic_carbocycles,numAliphaticCarbocycles,u32;
 num_saturated_heterocycles,numSaturatedHeterocycles,u32;num_saturated_carbocycles,numSaturatedCarbocycles,u32;
 num_amide_bonds,numAmideBonds,u32;num_atom_stereo_centers,numAtomStereoCenters,u32;
 num_unspecified_atom_stereo_centers,numUnspecifiedAtomStereoCenters,u32;
 num_spiro_atoms,numSpiroAtoms,u32;num_bridgehead_atoms,numBridgeheadAtoms,u32;qed,qed,f64;
 max_partial_charge,maxPartialCharge,f64;min_partial_charge,minPartialCharge,f64;
 max_abs_partial_charge,maxAbsPartialCharge,f64;min_abs_partial_charge,minAbsPartialCharge,f64;
 slogp_vsa_1,slogpVsa1,f64;slogp_vsa_2,slogpVsa2,f64;slogp_vsa_3,slogpVsa3,f64;
 slogp_vsa_4,slogpVsa4,f64;slogp_vsa_5,slogpVsa5,f64;slogp_vsa_6,slogpVsa6,f64;
 slogp_vsa_7,slogpVsa7,f64;slogp_vsa_8,slogpVsa8,f64;slogp_vsa_9,slogpVsa9,f64;
 slogp_vsa_10,slogpVsa10,f64;slogp_vsa_11,slogpVsa11,f64;slogp_vsa_12,slogpVsa12,f64;
 smr_vsa_1,smrVsa1,f64;smr_vsa_2,smrVsa2,f64;smr_vsa_3,smrVsa3,f64;
 smr_vsa_4,smrVsa4,f64;smr_vsa_5,smrVsa5,f64;smr_vsa_6,smrVsa6,f64;
 smr_vsa_7,smrVsa7,f64;smr_vsa_8,smrVsa8,f64;smr_vsa_9,smrVsa9,f64;smr_vsa_10,smrVsa10,f64;
 peoe_vsa_1,peoeVsa1,f64;peoe_vsa_2,peoeVsa2,f64;peoe_vsa_3,peoeVsa3,f64;
 peoe_vsa_4,peoeVsa4,f64;peoe_vsa_5,peoeVsa5,f64;peoe_vsa_6,peoeVsa6,f64;
 peoe_vsa_7,peoeVsa7,f64;peoe_vsa_8,peoeVsa8,f64;peoe_vsa_9,peoeVsa9,f64;
 peoe_vsa_10,peoeVsa10,f64;peoe_vsa_11,peoeVsa11,f64;peoe_vsa_12,peoeVsa12,f64;
 peoe_vsa_13,peoeVsa13,f64;peoe_vsa_14,peoeVsa14,f64;
}
macro_rules! force_queries{($($name:ident,$native:ident,$js:ident);+$(;)?)=>{paste::paste!{#[wasm_bindgen]impl Molecule{$(
    #[wasm_bindgen(js_name=$js)]pub fn $name(&self,force:Option<tsify::Ts<input::BoolInput>>)->Result<f64,JsValue>{self.0.$native(force.map(input::read::<bool>).transpose()?.unwrap_or(false)).map_err(DescriptorReadError::from_rust)}
    #[wasm_bindgen(js_name=[<$js WithParams>])]pub fn $native(&self,force:tsify::Ts<input::BoolInput>)->Result<f64,JsValue>{self.$name(Some(force))}
)+}}};}
force_queries! {
 chi_0_v,chi_0_v_with_params,chi0V;chi_1_v,chi_1_v_with_params,chi1V;chi_2_v,chi_2_v_with_params,chi2V;chi_3_v,chi_3_v_with_params,chi3V;chi_4_v,chi_4_v_with_params,chi4V;
 chi_0_n,chi_0_n_with_params,chi0N;chi_1_n,chi_1_n_with_params,chi1N;chi_2_n,chi_2_n_with_params,chi2N;chi_3_n,chi_3_n_with_params,chi3N;chi_4_n,chi_4_n_with_params,chi4N;
}
#[wasm_bindgen]
impl Molecule {
    pub fn mqns(&self, force: Option<tsify::Ts<input::BoolInput>>) -> Result<Vec<u32>, JsValue> {
        self.0
            .mqns(force.map(input::read::<bool>).transpose()?.unwrap_or(false))
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=chiNV)]
    pub fn chi_n_v(
        &self,
        order: tsify::Ts<input::U32Input>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .chi_n_v_with_params(
                input::read(order)?,
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=chiNN)]
    pub fn chi_n_n(
        &self,
        order: tsify::Ts<input::U32Input>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .chi_n_n_with_params(
                input::read(order)?,
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=molecularWeight)]
    pub fn molecular_weight(
        &self,
        only_heavy: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .molecular_weight_with_params(
                only_heavy
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=exactMolecularWeight)]
    pub fn exact_molecular_weight(
        &self,
        only_heavy: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .exact_molecular_weight_with_params(
                only_heavy
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=molecularFormula)]
    pub fn molecular_formula(
        &self,
        separate_isotopes: Option<tsify::Ts<input::BoolInput>>,
        abbreviate_h_isotopes: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<String, JsValue> {
        self.0
            .molecular_formula_with_params(
                separate_isotopes
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(false),
                abbreviate_h_isotopes
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    pub fn tpsa(
        &self,
        include_sulfur_phosphorus: Option<tsify::Ts<input::BoolInput>>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .tpsa_with_params(
                include_sulfur_phosphorus
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(false),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=labuteAsa)]
    pub fn labute_asa(
        &self,
        include_hydrogens: Option<tsify::Ts<input::BoolInput>>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<f64, JsValue> {
        self.0
            .labute_asa_with_params(
                include_hydrogens
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(true),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=slogpVsa)]
    pub fn slogp_vsa(
        &self,
        bins: Option<<Vec<f64> as Field>::Input>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<Vec<f64>, JsValue> {
        let bins = bins.map(Vec::<f64>::input).transpose()?;
        self.0
            .slogp_vsa_with_params(
                bins.as_deref(),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
    #[wasm_bindgen(js_name=smrVsa)]
    pub fn smr_vsa(
        &self,
        bins: Option<<Vec<f64> as Field>::Input>,
        force: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<Vec<f64>, JsValue> {
        let bins = bins.map(Vec::<f64>::input).transpose()?;
        self.0
            .smr_vsa_with_params(
                bins.as_deref(),
                force.map(input::read::<bool>).transpose()?.unwrap_or(false),
            )
            .map_err(DescriptorReadError::from_rust)
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=numRotatableBondsWithParams)]
    pub fn num_rotatable_bonds_with_params(
        &self,
        options: crate::enums::RotatableBondsOptionsInput,
    ) -> Result<u32, JsValue> {
        self.num_rotatable_bonds(Some(options))
    }
    #[wasm_bindgen(js_name=crippenDescriptorsWithParams)]
    pub fn crippen_descriptors_with_params(
        &self,
        include_hydrogens: tsify::Ts<input::BoolInput>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<CrippenTotals, JsValue> {
        self.crippen_descriptors(Some(include_hydrogens), Some(force))
    }
    #[wasm_bindgen(js_name=labuteAsaContributionsWithParams)]
    pub fn labute_asa_contributions_with_params(
        &self,
        include_hydrogens: tsify::Ts<input::BoolInput>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<LabuteAsaContributions, JsValue> {
        self.labute_asa_contributions(Some(include_hydrogens), Some(force))
    }
    #[wasm_bindgen(js_name=chiNVWithParams)]
    pub fn chi_n_v_with_params(
        &self,
        order: tsify::Ts<input::U32Input>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.chi_n_v(order, Some(force))
    }
    #[wasm_bindgen(js_name=chiNNWithParams)]
    pub fn chi_n_n_with_params(
        &self,
        order: tsify::Ts<input::U32Input>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.chi_n_n(order, Some(force))
    }
    #[wasm_bindgen(js_name=molecularWeightWithParams)]
    pub fn molecular_weight_with_params(
        &self,
        only_heavy: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.molecular_weight(Some(only_heavy))
    }
    #[wasm_bindgen(js_name=exactMolecularWeightWithParams)]
    pub fn exact_molecular_weight_with_params(
        &self,
        only_heavy: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.exact_molecular_weight(Some(only_heavy))
    }
    #[wasm_bindgen(js_name=molecularFormulaWithParams)]
    pub fn molecular_formula_with_params(
        &self,
        separate_isotopes: tsify::Ts<input::BoolInput>,
        abbreviate_h_isotopes: tsify::Ts<input::BoolInput>,
    ) -> Result<String, JsValue> {
        self.molecular_formula(Some(separate_isotopes), Some(abbreviate_h_isotopes))
    }
    #[wasm_bindgen(js_name=tpsaWithParams)]
    pub fn tpsa_with_params(
        &self,
        include_sulfur_phosphorus: tsify::Ts<input::BoolInput>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.tpsa(Some(include_sulfur_phosphorus), Some(force))
    }
    #[wasm_bindgen(js_name=labuteAsaWithParams)]
    pub fn labute_asa_with_params(
        &self,
        include_hydrogens: tsify::Ts<input::BoolInput>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<f64, JsValue> {
        self.labute_asa(Some(include_hydrogens), Some(force))
    }
    #[wasm_bindgen(js_name=slogpVsaWithParams)]
    pub fn slogp_vsa_with_params(
        &self,
        bins: Option<<Vec<f64> as Field>::Input>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<Vec<f64>, JsValue> {
        self.slogp_vsa(bins, Some(force))
    }
    #[wasm_bindgen(js_name=smrVsaWithParams)]
    pub fn smr_vsa_with_params(
        &self,
        bins: Option<<Vec<f64> as Field>::Input>,
        force: tsify::Ts<input::BoolInput>,
    ) -> Result<Vec<f64>, JsValue> {
        self.smr_vsa(bins, Some(force))
    }
}
