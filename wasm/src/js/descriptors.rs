//! Every registered descriptor delegates to the canonical public molecule.
use crate::{
    Molecule,
    descriptor_errors::descriptor_read_error,
    descriptor_values::{CrippenTotals, LabuteAsaContributions, rotatable_option},
    host_values::{bool_value, sequence, type_error, u32_value},
};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
fn parse_bins(v: &JsValue) -> Result<Option<Vec<f64>>, JsValue> {
    if v.is_null() {
        return Ok(None);
    }
    sequence(v, "bins")?
        .iter()
        .map(|v| v.as_f64().ok_or_else(|| type_error("bin")))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=chi0)]

    pub fn chi_0(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_0()
        self.inner
            .chi_0()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi1)]

    pub fn chi_1(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_1()
        self.inner
            .chi_1()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=hallKierAlpha)]

    pub fn hall_kier_alpha(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .hall_kier_alpha()
        self.inner
            .hall_kier_alpha()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=hallKierAlphaWithContributions)]
    #[wasm_bindgen(unchecked_return_type = "[number, Float64Array]")]
    pub fn hall_kier_alpha_with_contributions(&self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: .hall_kier_alpha_with_contributions()
        self.inner
            .hall_kier_alpha_with_contributions()
            .map(|(total, rows)| {
                let out = js_sys::Array::new();
                out.push(&total.into());
                out.push(&js_sys::Float64Array::from(rows.as_slice()));
                out.into()
            })
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=kappa1)]

    pub fn kappa_1(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .kappa_1()
        self.inner
            .kappa_1()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=kappa2)]

    pub fn kappa_2(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .kappa_2()
        self.inner
            .kappa_2()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=kappa3)]

    pub fn kappa_3(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .kappa_3()
        self.inner
            .kappa_3()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=phi)]

    pub fn phi(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .phi()
        self.inner
            .phi()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=mqns)]

    pub fn mqns(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] force: JsValue,
    ) -> Result<Vec<u32>, JsValue> {
        // COSMolKit❗✔️: .mqns(force)
        let force = if force.is_undefined() {
            false
        } else {
            bool_value(&force, "force")?
        };
        self.inner
            .mqns(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi0V)]

    pub fn chi_0_v(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_0_v()
        self.inner
            .chi_0_v()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi1V)]

    pub fn chi_1_v(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_1_v()
        self.inner
            .chi_1_v()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi2V)]

    pub fn chi_2_v(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_2_v()
        self.inner
            .chi_2_v()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi3V)]

    pub fn chi_3_v(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_3_v()
        self.inner
            .chi_3_v()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi4V)]

    pub fn chi_4_v(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_4_v()
        self.inner
            .chi_4_v()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chiNV)]

    pub fn chi_n_v(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] order: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_n_v(order)
        self.inner
            .chi_n_v(u32_value(&order, "order")?)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi0N)]

    pub fn chi_0_n(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_0_n()
        self.inner
            .chi_0_n()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi1N)]

    pub fn chi_1_n(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_1_n()
        self.inner
            .chi_1_n()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi2N)]

    pub fn chi_2_n(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_2_n()
        self.inner
            .chi_2_n()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi3N)]

    pub fn chi_3_n(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_3_n()
        self.inner
            .chi_3_n()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi4N)]

    pub fn chi_4_n(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_4_n()
        self.inner
            .chi_4_n()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chiNN)]

    pub fn chi_n_n(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] order: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_n_n(order)
        self.inner
            .chi_n_n(u32_value(&order, "order")?)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=molecularWeight)]

    pub fn molecular_weight(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .molecular_weight()
        self.inner
            .molecular_weight()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=exactMolecularWeight)]

    pub fn exact_molecular_weight(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .exact_molecular_weight()
        self.inner
            .exact_molecular_weight()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=molecularFormula)]

    pub fn molecular_formula(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: .molecular_formula()
        self.inner
            .molecular_formula()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAmideBonds)]

    pub fn num_amide_bonds(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_amide_bonds()
        self.inner
            .num_amide_bonds()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numSpiroAtoms)]

    pub fn num_spiro_atoms(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_spiro_atoms()
        self.inner
            .num_spiro_atoms()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numBridgeheadAtoms)]

    pub fn num_bridgehead_atoms(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_bridgehead_atoms()
        self.inner
            .num_bridgehead_atoms()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAtomStereoCenters)]

    pub fn num_atom_stereo_centers(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_atom_stereo_centers()
        self.inner
            .num_atom_stereo_centers()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numUnspecifiedAtomStereoCenters)]

    pub fn num_unspecified_atom_stereo_centers(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_unspecified_atom_stereo_centers()
        self.inner
            .num_unspecified_atom_stereo_centers()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numRotatableBonds)]

    pub fn num_rotatable_bonds(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_rotatable_bonds()
        self.inner
            .num_rotatable_bonds()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numRotatableBondsWithParams)]

    pub fn num_rotatable_bonds_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "RotatableBondsOptions")] params: JsValue,
    ) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_rotatable_bonds_with_params(params)
        self.inner
            .num_rotatable_bonds_with_params(&rotatable_option(&params)?)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=molecularWeightWithParams)]

    pub fn molecular_weight_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] only_heavy: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .molecular_weight_with_params(only_heavy)
        let only_heavy = bool_value(&only_heavy, "only_heavy")?;
        self.inner
            .molecular_weight_with_params(only_heavy)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=exactMolecularWeightWithParams)]

    pub fn exact_molecular_weight_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] only_heavy: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .exact_molecular_weight_with_params(only_heavy)
        let only_heavy = bool_value(&only_heavy, "only_heavy")?;
        self.inner
            .exact_molecular_weight_with_params(only_heavy)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=molecularFormulaWithParams)]

    pub fn molecular_formula_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] separate_isotopes: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] abbreviate_h_isotopes: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
        let separate_isotopes = bool_value(&separate_isotopes, "separate_isotopes")?;
        let abbreviate_h_isotopes = bool_value(&abbreviate_h_isotopes, "abbreviate_h_isotopes")?;
        self.inner
            .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numHeavyAtoms)]

    pub fn num_heavy_atoms(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_heavy_atoms()
        self.inner
            .num_heavy_atoms()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=totalAtomCount)]

    pub fn total_atom_count(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .total_atom_count()
        self.inner
            .total_atom_count()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numRings)]

    pub fn num_rings(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_rings()
        self.inner
            .num_rings()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numHeterocycles)]

    pub fn num_heterocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_heterocycles()
        self.inner
            .num_heterocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numHeteroatoms)]

    pub fn num_heteroatoms(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_heteroatoms()
        self.inner
            .num_heteroatoms()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numHba)]

    pub fn num_hba(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_hba()
        self.inner
            .num_hba()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numHbd)]

    pub fn num_hbd(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_hbd()
        self.inner
            .num_hbd()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAromaticRings)]

    pub fn num_aromatic_rings(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aromatic_rings()
        self.inner
            .num_aromatic_rings()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numSaturatedRings)]

    pub fn num_saturated_rings(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_saturated_rings()
        self.inner
            .num_saturated_rings()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAliphaticRings)]

    pub fn num_aliphatic_rings(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aliphatic_rings()
        self.inner
            .num_aliphatic_rings()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAromaticHeterocycles)]

    pub fn num_aromatic_heterocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aromatic_heterocycles()
        self.inner
            .num_aromatic_heterocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAromaticCarbocycles)]

    pub fn num_aromatic_carbocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aromatic_carbocycles()
        self.inner
            .num_aromatic_carbocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAliphaticHeterocycles)]

    pub fn num_aliphatic_heterocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aliphatic_heterocycles()
        self.inner
            .num_aliphatic_heterocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAliphaticCarbocycles)]

    pub fn num_aliphatic_carbocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_aliphatic_carbocycles()
        self.inner
            .num_aliphatic_carbocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numSaturatedHeterocycles)]

    pub fn num_saturated_heterocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_saturated_heterocycles()
        self.inner
            .num_saturated_heterocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numSaturatedCarbocycles)]

    pub fn num_saturated_carbocycles(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .num_saturated_carbocycles()
        self.inner
            .num_saturated_carbocycles()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=lipinskiHba)]

    pub fn lipinski_hba(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .lipinski_hba()
        self.inner
            .lipinski_hba()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=lipinskiHbd)]

    pub fn lipinski_hbd(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: .lipinski_hbd()
        self.inner
            .lipinski_hbd()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fractionCsp3)]

    pub fn fraction_csp3(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .fraction_csp3()
        self.inner
            .fraction_csp3()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=crippenDescriptors)]

    pub fn crippen_descriptors(&self) -> Result<CrippenTotals, JsValue> {
        // COSMolKit❗✔️: .crippen_descriptors()
        self.inner
            .crippen_descriptors()
            .map(|inner| CrippenTotals { inner })
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=labuteAsa)]

    pub fn labute_asa(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .labute_asa()
        self.inner
            .labute_asa()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=labuteAsaContributions)]

    pub fn labute_asa_contributions(&self) -> Result<LabuteAsaContributions, JsValue> {
        // COSMolKit❗✔️: .labute_asa_contributions()
        self.inner
            .labute_asa_contributions()
            .map(|inner| LabuteAsaContributions { inner })
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=tpsa)]

    pub fn tpsa(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .tpsa()
        self.inner
            .tpsa()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa)]

    pub fn slogp_vsa(&self) -> Result<Vec<f64>, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa()
        self.inner
            .slogp_vsa()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa)]

    pub fn smr_vsa(&self) -> Result<Vec<f64>, JsValue> {
        // COSMolKit❗✔️: .smr_vsa()
        self.inner
            .smr_vsa()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa1)]

    pub fn slogp_vsa_1(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_1()
        self.inner
            .slogp_vsa_1()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa2)]

    pub fn slogp_vsa_2(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_2()
        self.inner
            .slogp_vsa_2()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa3)]

    pub fn slogp_vsa_3(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_3()
        self.inner
            .slogp_vsa_3()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa4)]

    pub fn slogp_vsa_4(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_4()
        self.inner
            .slogp_vsa_4()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa5)]

    pub fn slogp_vsa_5(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_5()
        self.inner
            .slogp_vsa_5()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa6)]

    pub fn slogp_vsa_6(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_6()
        self.inner
            .slogp_vsa_6()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa7)]

    pub fn slogp_vsa_7(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_7()
        self.inner
            .slogp_vsa_7()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa8)]

    pub fn slogp_vsa_8(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_8()
        self.inner
            .slogp_vsa_8()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa9)]

    pub fn slogp_vsa_9(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_9()
        self.inner
            .slogp_vsa_9()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa10)]

    pub fn slogp_vsa_10(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_10()
        self.inner
            .slogp_vsa_10()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa11)]

    pub fn slogp_vsa_11(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_11()
        self.inner
            .slogp_vsa_11()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsa12)]

    pub fn slogp_vsa_12(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_12()
        self.inner
            .slogp_vsa_12()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa1)]

    pub fn smr_vsa_1(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_1()
        self.inner
            .smr_vsa_1()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa2)]

    pub fn smr_vsa_2(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_2()
        self.inner
            .smr_vsa_2()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa3)]

    pub fn smr_vsa_3(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_3()
        self.inner
            .smr_vsa_3()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa4)]

    pub fn smr_vsa_4(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_4()
        self.inner
            .smr_vsa_4()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa5)]

    pub fn smr_vsa_5(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_5()
        self.inner
            .smr_vsa_5()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa6)]

    pub fn smr_vsa_6(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_6()
        self.inner
            .smr_vsa_6()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa7)]

    pub fn smr_vsa_7(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_7()
        self.inner
            .smr_vsa_7()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa8)]

    pub fn smr_vsa_8(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_8()
        self.inner
            .smr_vsa_8()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa9)]

    pub fn smr_vsa_9(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_9()
        self.inner
            .smr_vsa_9()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsa10)]

    pub fn smr_vsa_10(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_10()
        self.inner
            .smr_vsa_10()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=crippenDescriptorsWithParams)]

    pub fn crippen_descriptors_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] include_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<CrippenTotals, JsValue> {
        // COSMolKit❗✔️: .crippen_descriptors_with_params(include_hydrogens, force)
        let include_hydrogens = bool_value(&include_hydrogens, "include_hydrogens")?;
        let force = bool_value(&force, "force")?;
        self.inner
            .crippen_descriptors_with_params(include_hydrogens, force)
            .map(|inner| CrippenTotals { inner })
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=labuteAsaWithParams)]

    pub fn labute_asa_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] include_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .labute_asa_with_params(include_hydrogens, force)
        let include_hydrogens = bool_value(&include_hydrogens, "include_hydrogens")?;
        let force = bool_value(&force, "force")?;
        self.inner
            .labute_asa_with_params(include_hydrogens, force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=labuteAsaContributionsWithParams)]

    pub fn labute_asa_contributions_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] include_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<LabuteAsaContributions, JsValue> {
        // COSMolKit❗✔️: .labute_asa_contributions_with_params(include_hydrogens, force)
        let include_hydrogens = bool_value(&include_hydrogens, "include_hydrogens")?;
        let force = bool_value(&force, "force")?;
        self.inner
            .labute_asa_contributions_with_params(include_hydrogens, force)
            .map(|inner| LabuteAsaContributions { inner })
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=tpsaWithParams)]

    pub fn tpsa_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] include_sulfur_phosphorus: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .tpsa_with_params(include_sulfur_phosphorus, force)
        let include_sulfur_phosphorus =
            bool_value(&include_sulfur_phosphorus, "include_sulfur_phosphorus")?;
        let force = bool_value(&force, "force")?;
        self.inner
            .tpsa_with_params(include_sulfur_phosphorus, force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=slogpVsaWithParams)]

    pub fn slogp_vsa_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] bins: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<Vec<f64>, JsValue> {
        // COSMolKit❗✔️: .slogp_vsa_with_params(bins, force)
        let bins = parse_bins(&bins)?;
        let force = bool_value(&force, "force")?;
        self.inner
            .slogp_vsa_with_params(bins.as_deref(), force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=smrVsaWithParams)]

    pub fn smr_vsa_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] bins: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<Vec<f64>, JsValue> {
        // COSMolKit❗✔️: .smr_vsa_with_params(bins, force)
        let bins = parse_bins(&bins)?;
        let force = bool_value(&force, "force")?;
        self.inner
            .smr_vsa_with_params(bins.as_deref(), force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=qed)]

    pub fn qed(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .qed()
        self.inner
            .qed()
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi0VWithParams)]

    pub fn chi_0_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_0_v_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_0_v_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi1VWithParams)]

    pub fn chi_1_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_1_v_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_1_v_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi2VWithParams)]

    pub fn chi_2_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_2_v_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_2_v_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi3VWithParams)]

    pub fn chi_3_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_3_v_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_3_v_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi4VWithParams)]

    pub fn chi_4_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_4_v_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_4_v_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chiNVWithParams)]

    pub fn chi_n_v_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] order: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_n_v_with_params(order, force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_n_v_with_params(u32_value(&order, "order")?, force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi0NWithParams)]

    pub fn chi_0_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_0_n_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_0_n_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi1NWithParams)]

    pub fn chi_1_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_1_n_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_1_n_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi2NWithParams)]

    pub fn chi_2_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_2_n_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_2_n_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi3NWithParams)]

    pub fn chi_3_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_3_n_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_3_n_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chi4NWithParams)]

    pub fn chi_4_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_4_n_with_params(force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_4_n_with_params(force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=chiNNWithParams)]

    pub fn chi_n_n_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] order: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] force: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .chi_n_n_with_params(order, force)
        let force = bool_value(&force, "force")?;
        self.inner
            .chi_n_n_with_params(u32_value(&order, "order")?, force)
            .map_err(|e| descriptor_read_error(&e).unwrap_or_else(|e| e))
    }
}
