//! One public facade owner for all descriptor queries; no binding chemistry.
use super::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn chi_0(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_0()
        self.inner.borrow().chi_0()
    }
    pub fn chi_1(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_1()
        self.inner.borrow().chi_1()
    }
    pub fn hall_kier_alpha(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .hall_kier_alpha()
        self.inner.borrow().hall_kier_alpha()
    }
    pub fn hall_kier_alpha_with_contributions(
        &self,
    ) -> Result<(f64, Vec<f64>), ck::DescriptorReadError> {
        // COSMolKit❗✔️: .hall_kier_alpha_with_contributions()
        self.inner.borrow().hall_kier_alpha_with_contributions()
    }
    pub fn kappa_1(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .kappa_1()
        self.inner.borrow().kappa_1()
    }
    pub fn kappa_2(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .kappa_2()
        self.inner.borrow().kappa_2()
    }
    pub fn kappa_3(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .kappa_3()
        self.inner.borrow().kappa_3()
    }
    pub fn phi(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .phi()
        self.inner.borrow().phi()
    }
    pub fn mqns(&self, force: bool) -> Result<Vec<u32>, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .mqns(force)
        self.inner.borrow().mqns(force)
    }
    pub fn chi_0_v(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_0_v()
        self.inner.borrow().chi_0_v()
    }
    pub fn chi_1_v(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_1_v()
        self.inner.borrow().chi_1_v()
    }
    pub fn chi_2_v(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_2_v()
        self.inner.borrow().chi_2_v()
    }
    pub fn chi_3_v(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_3_v()
        self.inner.borrow().chi_3_v()
    }
    pub fn chi_4_v(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_4_v()
        self.inner.borrow().chi_4_v()
    }
    pub fn chi_n_v(&self, order: u32) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_n_v(order)
        self.inner.borrow().chi_n_v(order)
    }
    pub fn chi_0_n(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_0_n()
        self.inner.borrow().chi_0_n()
    }
    pub fn chi_1_n(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_1_n()
        self.inner.borrow().chi_1_n()
    }
    pub fn chi_2_n(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_2_n()
        self.inner.borrow().chi_2_n()
    }
    pub fn chi_3_n(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_3_n()
        self.inner.borrow().chi_3_n()
    }
    pub fn chi_4_n(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_4_n()
        self.inner.borrow().chi_4_n()
    }
    pub fn chi_n_n(&self, order: u32) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_n_n(order)
        self.inner.borrow().chi_n_n(order)
    }
    pub fn molecular_weight(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .molecular_weight()
        self.inner.borrow().molecular_weight()
    }
    pub fn exact_molecular_weight(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .exact_molecular_weight()
        self.inner.borrow().exact_molecular_weight()
    }
    pub fn molecular_formula(&self) -> Result<String, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .molecular_formula()
        self.inner.borrow().molecular_formula()
    }
    pub fn num_amide_bonds(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_amide_bonds()
        self.inner.borrow().num_amide_bonds()
    }
    pub fn num_spiro_atoms(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_spiro_atoms()
        self.inner.borrow().num_spiro_atoms()
    }
    pub fn num_bridgehead_atoms(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_bridgehead_atoms()
        self.inner.borrow().num_bridgehead_atoms()
    }
    pub fn num_atom_stereo_centers(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_atom_stereo_centers()
        self.inner.borrow().num_atom_stereo_centers()
    }
    pub fn num_unspecified_atom_stereo_centers(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_unspecified_atom_stereo_centers()
        self.inner.borrow().num_unspecified_atom_stereo_centers()
    }
    pub fn num_rotatable_bonds(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_rotatable_bonds()
        self.inner.borrow().num_rotatable_bonds()
    }
    pub fn num_rotatable_bonds_with_params(
        &self,
        params: &ck::RotatableBondsOptions,
    ) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_rotatable_bonds_with_params(params)
        self.inner.borrow().num_rotatable_bonds_with_params(params)
    }
    pub fn molecular_weight_with_params(
        &self,
        only_heavy: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .molecular_weight_with_params(only_heavy)
        self.inner.borrow().molecular_weight_with_params(only_heavy)
    }
    pub fn exact_molecular_weight_with_params(
        &self,
        only_heavy: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .exact_molecular_weight_with_params(only_heavy)
        self.inner
            .borrow()
            .exact_molecular_weight_with_params(only_heavy)
    }
    pub fn molecular_formula_with_params(
        &self,
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    ) -> Result<String, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
        self.inner
            .borrow()
            .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
    }
    pub fn num_heavy_atoms(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_heavy_atoms()
        self.inner.borrow().num_heavy_atoms()
    }
    pub fn total_atom_count(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .total_atom_count()
        self.inner.borrow().total_atom_count()
    }
    pub fn num_rings(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_rings()
        self.inner.borrow().num_rings()
    }
    pub fn num_heterocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_heterocycles()
        self.inner.borrow().num_heterocycles()
    }
    pub fn num_heteroatoms(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_heteroatoms()
        self.inner.borrow().num_heteroatoms()
    }
    pub fn num_hba(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_hba()
        self.inner.borrow().num_hba()
    }
    pub fn num_hbd(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_hbd()
        self.inner.borrow().num_hbd()
    }
    pub fn num_aromatic_rings(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aromatic_rings()
        self.inner.borrow().num_aromatic_rings()
    }
    pub fn num_saturated_rings(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_saturated_rings()
        self.inner.borrow().num_saturated_rings()
    }
    pub fn num_aliphatic_rings(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aliphatic_rings()
        self.inner.borrow().num_aliphatic_rings()
    }
    pub fn num_aromatic_heterocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aromatic_heterocycles()
        self.inner.borrow().num_aromatic_heterocycles()
    }
    pub fn num_aromatic_carbocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aromatic_carbocycles()
        self.inner.borrow().num_aromatic_carbocycles()
    }
    pub fn num_aliphatic_heterocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aliphatic_heterocycles()
        self.inner.borrow().num_aliphatic_heterocycles()
    }
    pub fn num_aliphatic_carbocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_aliphatic_carbocycles()
        self.inner.borrow().num_aliphatic_carbocycles()
    }
    pub fn num_saturated_heterocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_saturated_heterocycles()
        self.inner.borrow().num_saturated_heterocycles()
    }
    pub fn num_saturated_carbocycles(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .num_saturated_carbocycles()
        self.inner.borrow().num_saturated_carbocycles()
    }
    pub fn lipinski_hba(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .lipinski_hba()
        self.inner.borrow().lipinski_hba()
    }
    pub fn lipinski_hbd(&self) -> Result<u32, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .lipinski_hbd()
        self.inner.borrow().lipinski_hbd()
    }
    pub fn fraction_csp3(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .fraction_csp3()
        self.inner.borrow().fraction_csp3()
    }
    pub fn crippen_descriptors(&self) -> Result<ck::CrippenTotals, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .crippen_descriptors()
        self.inner.borrow().crippen_descriptors()
    }
    pub fn labute_asa(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .labute_asa()
        self.inner.borrow().labute_asa()
    }
    pub fn labute_asa_contributions(
        &self,
    ) -> Result<ck::LabuteAsaContributions, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .labute_asa_contributions()
        self.inner.borrow().labute_asa_contributions()
    }
    pub fn tpsa(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .tpsa()
        self.inner.borrow().tpsa()
    }
    pub fn slogp_vsa(&self) -> Result<Vec<f64>, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa()
        self.inner.borrow().slogp_vsa()
    }
    pub fn smr_vsa(&self) -> Result<Vec<f64>, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa()
        self.inner.borrow().smr_vsa()
    }
    pub fn slogp_vsa_1(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_1()
        self.inner.borrow().slogp_vsa_1()
    }
    pub fn slogp_vsa_2(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_2()
        self.inner.borrow().slogp_vsa_2()
    }
    pub fn slogp_vsa_3(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_3()
        self.inner.borrow().slogp_vsa_3()
    }
    pub fn slogp_vsa_4(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_4()
        self.inner.borrow().slogp_vsa_4()
    }
    pub fn slogp_vsa_5(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_5()
        self.inner.borrow().slogp_vsa_5()
    }
    pub fn slogp_vsa_6(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_6()
        self.inner.borrow().slogp_vsa_6()
    }
    pub fn slogp_vsa_7(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_7()
        self.inner.borrow().slogp_vsa_7()
    }
    pub fn slogp_vsa_8(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_8()
        self.inner.borrow().slogp_vsa_8()
    }
    pub fn slogp_vsa_9(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_9()
        self.inner.borrow().slogp_vsa_9()
    }
    pub fn slogp_vsa_10(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_10()
        self.inner.borrow().slogp_vsa_10()
    }
    pub fn slogp_vsa_11(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_11()
        self.inner.borrow().slogp_vsa_11()
    }
    pub fn slogp_vsa_12(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_12()
        self.inner.borrow().slogp_vsa_12()
    }
    pub fn smr_vsa_1(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_1()
        self.inner.borrow().smr_vsa_1()
    }
    pub fn smr_vsa_2(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_2()
        self.inner.borrow().smr_vsa_2()
    }
    pub fn smr_vsa_3(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_3()
        self.inner.borrow().smr_vsa_3()
    }
    pub fn smr_vsa_4(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_4()
        self.inner.borrow().smr_vsa_4()
    }
    pub fn smr_vsa_5(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_5()
        self.inner.borrow().smr_vsa_5()
    }
    pub fn smr_vsa_6(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_6()
        self.inner.borrow().smr_vsa_6()
    }
    pub fn smr_vsa_7(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_7()
        self.inner.borrow().smr_vsa_7()
    }
    pub fn smr_vsa_8(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_8()
        self.inner.borrow().smr_vsa_8()
    }
    pub fn smr_vsa_9(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_9()
        self.inner.borrow().smr_vsa_9()
    }
    pub fn smr_vsa_10(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_10()
        self.inner.borrow().smr_vsa_10()
    }
    pub fn crippen_descriptors_with_params(
        &self,
        include_hydrogens: bool,
        force: bool,
    ) -> Result<ck::CrippenTotals, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .crippen_descriptors_with_params(include_hydrogens, force)
        self.inner
            .borrow()
            .crippen_descriptors_with_params(include_hydrogens, force)
    }
    pub fn labute_asa_with_params(
        &self,
        include_hydrogens: bool,
        force: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .labute_asa_with_params(include_hydrogens, force)
        self.inner
            .borrow()
            .labute_asa_with_params(include_hydrogens, force)
    }
    pub fn labute_asa_contributions_with_params(
        &self,
        include_hydrogens: bool,
        force: bool,
    ) -> Result<ck::LabuteAsaContributions, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .labute_asa_contributions_with_params(include_hydrogens, force)
        self.inner
            .borrow()
            .labute_asa_contributions_with_params(include_hydrogens, force)
    }
    pub fn tpsa_with_params(
        &self,
        include_sulfur_phosphorus: bool,
        force: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .tpsa_with_params(include_sulfur_phosphorus, force)
        self.inner
            .borrow()
            .tpsa_with_params(include_sulfur_phosphorus, force)
    }
    pub fn slogp_vsa_with_params(
        &self,
        bins: Option<&[f64]>,
        force: bool,
    ) -> Result<Vec<f64>, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .slogp_vsa_with_params(bins, force)
        self.inner.borrow().slogp_vsa_with_params(bins, force)
    }
    pub fn smr_vsa_with_params(
        &self,
        bins: Option<&[f64]>,
        force: bool,
    ) -> Result<Vec<f64>, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .smr_vsa_with_params(bins, force)
        self.inner.borrow().smr_vsa_with_params(bins, force)
    }
    pub fn qed(&self) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .qed()
        self.inner.borrow().qed()
    }
    pub fn chi_0_v_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_0_v_with_params(force)
        self.inner.borrow().chi_0_v_with_params(force)
    }
    pub fn chi_1_v_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_1_v_with_params(force)
        self.inner.borrow().chi_1_v_with_params(force)
    }
    pub fn chi_2_v_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_2_v_with_params(force)
        self.inner.borrow().chi_2_v_with_params(force)
    }
    pub fn chi_3_v_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_3_v_with_params(force)
        self.inner.borrow().chi_3_v_with_params(force)
    }
    pub fn chi_4_v_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_4_v_with_params(force)
        self.inner.borrow().chi_4_v_with_params(force)
    }
    pub fn chi_n_v_with_params(
        &self,
        order: u32,
        force: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_n_v_with_params(order, force)
        self.inner.borrow().chi_n_v_with_params(order, force)
    }
    pub fn chi_0_n_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_0_n_with_params(force)
        self.inner.borrow().chi_0_n_with_params(force)
    }
    pub fn chi_1_n_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_1_n_with_params(force)
        self.inner.borrow().chi_1_n_with_params(force)
    }
    pub fn chi_2_n_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_2_n_with_params(force)
        self.inner.borrow().chi_2_n_with_params(force)
    }
    pub fn chi_3_n_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_3_n_with_params(force)
        self.inner.borrow().chi_3_n_with_params(force)
    }
    pub fn chi_4_n_with_params(&self, force: bool) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_4_n_with_params(force)
        self.inner.borrow().chi_4_n_with_params(force)
    }
    pub fn chi_n_n_with_params(
        &self,
        order: u32,
        force: bool,
    ) -> Result<f64, ck::DescriptorReadError> {
        // COSMolKit❗✔️: .chi_n_n_with_params(order, force)
        self.inner.borrow().chi_n_n_with_params(order, force)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn descriptor_registry_all_calls_and_cache_semantics() {
        for text in ["CCO", "c1ccncc1", "F[C@](Cl)(Br)I"] {
            let m = Molecule::from_smiles(text).unwrap();
            let before = m.to_smiles().unwrap();
            {
                let v = m.chi_0().expect("chi_0");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_1().expect("chi_1");
                assert!(v.is_finite());
            }
            {
                let v = m.hall_kier_alpha().expect("hall_kier_alpha");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .hall_kier_alpha_with_contributions()
                    .expect("hall_kier_alpha_with_contributions");
                assert!(v.0.is_finite());
                assert_eq!(v.1.len(), usize::try_from(m.num_atoms()).unwrap());
            }
            {
                let v = m.kappa_1().expect("kappa_1");
                assert!(v.is_finite());
            }
            {
                let v = m.kappa_2().expect("kappa_2");
                assert!(v.is_finite());
            }
            {
                let v = m.kappa_3().expect("kappa_3");
                assert!(v.is_finite());
            }
            {
                let v = m.phi().expect("phi");
                assert!(v.is_finite());
            }
            {
                let v = m.mqns(false).expect("mqns");
                assert_eq!(v.len(), 42);
            }
            {
                let v = m.chi_0_v().expect("chi_0_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_1_v().expect("chi_1_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_2_v().expect("chi_2_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_3_v().expect("chi_3_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_4_v().expect("chi_4_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_n_v(2).expect("chi_n_v");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_0_n().expect("chi_0_n");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_1_n().expect("chi_1_n");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_2_n().expect("chi_2_n");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_3_n().expect("chi_3_n");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_4_n().expect("chi_4_n");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_n_n(2).expect("chi_n_n");
                assert!(v.is_finite());
            }
            {
                let v = m.molecular_weight().expect("molecular_weight");
                assert!(v.is_finite());
            }
            {
                let v = m.exact_molecular_weight().expect("exact_molecular_weight");
                assert!(v.is_finite());
            }
            {
                let v = m.molecular_formula().expect("molecular_formula");
                assert!(!v.is_empty());
            }
            {
                let v = m.num_amide_bonds().expect("num_amide_bonds");
                let _: u32 = v;
            }
            {
                let v = m.num_spiro_atoms().expect("num_spiro_atoms");
                let _: u32 = v;
            }
            {
                let v = m.num_bridgehead_atoms().expect("num_bridgehead_atoms");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_atom_stereo_centers()
                    .expect("num_atom_stereo_centers");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_unspecified_atom_stereo_centers()
                    .expect("num_unspecified_atom_stereo_centers");
                let _: u32 = v;
            }
            {
                let v = m.num_rotatable_bonds().expect("num_rotatable_bonds");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_rotatable_bonds_with_params(&ck::RotatableBondsOptions::Default)
                    .expect("num_rotatable_bonds_with_params");
                let _: u32 = v;
            }
            {
                let v = m
                    .molecular_weight_with_params(false)
                    .expect("molecular_weight_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .exact_molecular_weight_with_params(false)
                    .expect("exact_molecular_weight_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .molecular_formula_with_params(true, true)
                    .expect("molecular_formula_with_params");
                assert!(!v.is_empty());
            }
            {
                let v = m.num_heavy_atoms().expect("num_heavy_atoms");
                let _: u32 = v;
            }
            {
                let v = m.total_atom_count().expect("total_atom_count");
                let _: u32 = v;
            }
            {
                let v = m.num_rings().expect("num_rings");
                let _: u32 = v;
            }
            {
                let v = m.num_heterocycles().expect("num_heterocycles");
                let _: u32 = v;
            }
            {
                let v = m.num_heteroatoms().expect("num_heteroatoms");
                let _: u32 = v;
            }
            {
                let v = m.num_hba().expect("num_hba");
                let _: u32 = v;
            }
            {
                let v = m.num_hbd().expect("num_hbd");
                let _: u32 = v;
            }
            {
                let v = m.num_aromatic_rings().expect("num_aromatic_rings");
                let _: u32 = v;
            }
            {
                let v = m.num_saturated_rings().expect("num_saturated_rings");
                let _: u32 = v;
            }
            {
                let v = m.num_aliphatic_rings().expect("num_aliphatic_rings");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_aromatic_heterocycles()
                    .expect("num_aromatic_heterocycles");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_aromatic_carbocycles()
                    .expect("num_aromatic_carbocycles");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_aliphatic_heterocycles()
                    .expect("num_aliphatic_heterocycles");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_aliphatic_carbocycles()
                    .expect("num_aliphatic_carbocycles");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_saturated_heterocycles()
                    .expect("num_saturated_heterocycles");
                let _: u32 = v;
            }
            {
                let v = m
                    .num_saturated_carbocycles()
                    .expect("num_saturated_carbocycles");
                let _: u32 = v;
            }
            {
                let v = m.lipinski_hba().expect("lipinski_hba");
                let _: u32 = v;
            }
            {
                let v = m.lipinski_hbd().expect("lipinski_hbd");
                let _: u32 = v;
            }
            {
                let v = m.fraction_csp3().expect("fraction_csp3");
                assert!(v.is_finite());
            }
            {
                let v = m.crippen_descriptors().expect("crippen_descriptors");
                assert!(v.logp.is_finite());
                assert!(v.molar_refractivity.is_finite());
            }
            {
                let v = m.labute_asa().expect("labute_asa");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .labute_asa_contributions()
                    .expect("labute_asa_contributions");
                assert!(v.asa.is_finite());
                assert_eq!(
                    v.atom_contributions.len(),
                    usize::try_from(m.num_atoms()).unwrap()
                );
            }
            {
                let v = m.tpsa().expect("tpsa");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa().expect("slogp_vsa");
                assert!(v.iter().all(|x| x.is_finite()));
            }
            {
                let v = m.smr_vsa().expect("smr_vsa");
                assert!(v.iter().all(|x| x.is_finite()));
            }
            {
                let v = m.slogp_vsa_1().expect("slogp_vsa_1");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_2().expect("slogp_vsa_2");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_3().expect("slogp_vsa_3");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_4().expect("slogp_vsa_4");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_5().expect("slogp_vsa_5");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_6().expect("slogp_vsa_6");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_7().expect("slogp_vsa_7");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_8().expect("slogp_vsa_8");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_9().expect("slogp_vsa_9");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_10().expect("slogp_vsa_10");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_11().expect("slogp_vsa_11");
                assert!(v.is_finite());
            }
            {
                let v = m.slogp_vsa_12().expect("slogp_vsa_12");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_1().expect("smr_vsa_1");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_2().expect("smr_vsa_2");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_3().expect("smr_vsa_3");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_4().expect("smr_vsa_4");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_5().expect("smr_vsa_5");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_6().expect("smr_vsa_6");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_7().expect("smr_vsa_7");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_8().expect("smr_vsa_8");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_9().expect("smr_vsa_9");
                assert!(v.is_finite());
            }
            {
                let v = m.smr_vsa_10().expect("smr_vsa_10");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .crippen_descriptors_with_params(true, false)
                    .expect("crippen_descriptors_with_params");
                assert!(v.logp.is_finite());
                assert!(v.molar_refractivity.is_finite());
            }
            {
                let v = m
                    .labute_asa_with_params(true, false)
                    .expect("labute_asa_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .labute_asa_contributions_with_params(true, false)
                    .expect("labute_asa_contributions_with_params");
                assert!(v.asa.is_finite());
                assert_eq!(
                    v.atom_contributions.len(),
                    usize::try_from(m.num_atoms()).unwrap()
                );
            }
            {
                let v = m.tpsa_with_params(false, false).expect("tpsa_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .slogp_vsa_with_params(None, false)
                    .expect("slogp_vsa_with_params");
                assert!(v.iter().all(|x| x.is_finite()));
            }
            {
                let v = m
                    .smr_vsa_with_params(None, false)
                    .expect("smr_vsa_with_params");
                assert!(v.iter().all(|x| x.is_finite()));
            }
            {
                let v = m.qed().expect("qed");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_0_v_with_params(false).expect("chi_0_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_1_v_with_params(false).expect("chi_1_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_2_v_with_params(false).expect("chi_2_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_3_v_with_params(false).expect("chi_3_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_4_v_with_params(false).expect("chi_4_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .chi_n_v_with_params(2, false)
                    .expect("chi_n_v_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_0_n_with_params(false).expect("chi_0_n_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_1_n_with_params(false).expect("chi_1_n_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_2_n_with_params(false).expect("chi_2_n_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_3_n_with_params(false).expect("chi_3_n_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m.chi_4_n_with_params(false).expect("chi_4_n_with_params");
                assert!(v.is_finite());
            }
            {
                let v = m
                    .chi_n_n_with_params(2, false)
                    .expect("chi_n_n_with_params");
                assert!(v.is_finite());
            }
            assert_eq!(m.to_smiles().unwrap(), before);
            assert!(m.coordinates_2d().is_empty());
            assert_eq!(m.num_3d_conformers(), 0);
        }
        let m = Molecule::from_smiles("CCO").unwrap();
        let cold = m.crippen_descriptors_with_params(false, false).unwrap();
        assert_eq!(cold.logp.to_bits(), (-0.3487_f64).to_bits());
        assert_eq!(
            m.crippen_descriptors_with_params(true, false)
                .unwrap()
                .logp
                .to_bits(),
            cold.logp.to_bits()
        );
        let forced = m.crippen_descriptors_with_params(true, true).unwrap();
        assert_eq!(
            forced.logp.to_bits(),
            (-0.0014000000000000123_f64).to_bits()
        );
        assert_eq!(
            m.crippen_descriptors_with_params(false, false)
                .unwrap()
                .logp
                .to_bits(),
            forced.logp.to_bits()
        );
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(raw.molecular_weight().unwrap().is_finite());
        assert!(raw.exact_molecular_weight().unwrap().is_finite());
        assert_eq!(raw.molecular_formula().unwrap(), "C2H6O");
        assert!(matches!(
            raw.num_hba(),
            Err(ck::DescriptorReadError::MissingPreparedValence)
        ));
        assert!(matches!(
            raw.num_rings(),
            Err(ck::DescriptorReadError::MissingInitializedRings)
        ));
        assert_eq!(raw.num_aromatic_rings().unwrap(), 0);
        assert_eq!(raw.num_heteroatoms().unwrap(), 1);
    }
}
