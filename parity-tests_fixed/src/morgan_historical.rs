//! Historical ChEMBL Morgan calls through existing public generators.
use crate::molecular::{MorganDenseBitsOutput, observed_morgan_additional_output};
use cosmolkit::{
    FingerprintAdditionalOutput, Molecule, MorganAtomInvariantsGenerator,
    MorganBondInvariantsGenerator, MorganCallParams, MorganFingerprintGenerator, MorganParams,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Recipe {
    Default,
    Chiral,
    Radius3,
    NoBondTypes,
    CountSimulation,
    RedundantEnvironments,
    FeatureInvariants,
    ConnectivityWithoutRingMembership,
    ChiralBondInvariants,
    TwoBitsPerFeature,
    FromFirst,
    IgnoreFirst,
    CustomAtomInvariants,
    OnlyNonzeroCustomInvariants,
    CustomBondInvariants,
}

impl Recipe {
    pub const ALL: [Self; 15] = [
        Self::Default,
        Self::Chiral,
        Self::Radius3,
        Self::NoBondTypes,
        Self::CountSimulation,
        Self::RedundantEnvironments,
        Self::FeatureInvariants,
        Self::ConnectivityWithoutRingMembership,
        Self::ChiralBondInvariants,
        Self::TwoBitsPerFeature,
        Self::FromFirst,
        Self::IgnoreFirst,
        Self::CustomAtomInvariants,
        Self::OnlyNonzeroCustomInvariants,
        Self::CustomBondInvariants,
    ];

    pub fn eligible(self, atoms: usize, bonds: usize) -> bool {
        match self {
            Self::FromFirst
            | Self::IgnoreFirst
            | Self::CustomAtomInvariants
            | Self::OnlyNonzeroCustomInvariants => atoms > 0,
            Self::CustomBondInvariants => bonds > 0,
            _ => true,
        }
    }

    pub fn params(self) -> MorganParams {
        // The historical ChEMBL calls explicitly selected radius2. The
        // current source-shaped public generator default is radius3.
        let mut p = MorganParams {
            radius: 2,
            ..MorganParams::default()
        };
        match self {
            Self::Chiral => p.include_chirality = true,
            Self::Radius3 => p.radius = 3,
            Self::NoBondTypes => {
                p.fp_size = 1024;
                p.use_bond_types = false;
            }
            Self::CountSimulation => p.count_simulation = true,
            Self::RedundantEnvironments => p.include_redundant_environments = true,
            Self::TwoBitsPerFeature => p.bits_per_feature = 2,
            Self::OnlyNonzeroCustomInvariants => p.only_nonzero_invariants = true,
            _ => {}
        }
        p
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Observation {
    // Compare eligibility independently, so a wrong CK count cannot silently
    // remove a historical conditional recipe from the selected workload.
    pub atom_count: usize,
    pub bond_count: usize,
    pub fingerprint: MorganDenseBitsOutput,
}

pub fn run(mol: &Molecule, recipe: Recipe) -> Result<Observation, String> {
    if !recipe.eligible(mol.num_atoms(), mol.num_bonds()) {
        return Err("historical Morgan recipe does not meet its source eligibility".into());
    }
    let params = recipe.params();
    let atom_provider = match recipe {
        Recipe::FeatureInvariants => Some(MorganAtomInvariantsGenerator::features(None)),
        Recipe::ConnectivityWithoutRingMembership => {
            Some(MorganAtomInvariantsGenerator::connectivity(false))
        }
        _ => None,
    };
    let bond_provider = (recipe == Recipe::ChiralBondInvariants)
        .then(|| MorganBondInvariantsGenerator::new(true, true));
    let generator = MorganFingerprintGenerator::new(
        Some(&params),
        atom_provider.as_ref(),
        bond_provider.as_ref(),
    )
    .map_err(|e| e.to_string())?;
    let mut call = MorganCallParams::default();
    match recipe {
        Recipe::FromFirst => call.from_atoms = Some(vec![0]),
        Recipe::IgnoreFirst => call.ignore_atoms = Some(vec![0]),
        Recipe::CustomAtomInvariants => {
            call.custom_atom_invariants =
                Some((0..mol.num_atoms() as u32).map(|i| i + 1).collect());
        }
        Recipe::OnlyNonzeroCustomInvariants => {
            call.custom_atom_invariants = Some(
                (0..mol.num_atoms() as u32)
                    .map(|i| if i % 2 == 0 { 0 } else { i + 1 })
                    .collect(),
            );
        }
        Recipe::CustomBondInvariants => {
            call.custom_bond_invariants =
                Some((0..mol.num_bonds() as u32).map(|i| i + 7).collect());
        }
        _ => {}
    }
    let mut extra = FingerprintAdditionalOutput::default();
    extra.allocate_atom_counts();
    extra.allocate_atom_to_bits();
    extra.allocate_bit_info_map();
    extra.allocate_bit_paths();
    extra.allocate_atoms_per_bit();
    let fp = mol
        .fingerprint_morgan_with_generator(&generator, Some(&call), Some(&mut extra))
        .map_err(|e| e.to_string())?;
    Ok(Observation {
        atom_count: mol.num_atoms(),
        bond_count: mol.num_bonds(),
        fingerprint: MorganDenseBitsOutput {
            length: fp.n_bits(),
            on_bits: fp.on_bits(),
            additional_output: observed_morgan_additional_output(&extra),
        },
    })
}

pub fn valid(recipe: Recipe, value: &Observation) -> bool {
    let fp = &value.fingerprint;
    recipe.eligible(value.atom_count, value.bond_count)
        && fp.length == recipe.params().fp_size
        && fp.on_bits.iter().all(|&bit| bit < fp.length)
        && fp.on_bits.windows(2).all(|pair| pair[0] < pair[1])
        && crate::molecular::valid_morgan_additional_output(&fp.additional_output)
        && fp
            .additional_output
            .atom_counts
            .as_ref()
            .is_some_and(|v| v.len() == value.atom_count)
        && fp
            .additional_output
            .atom_to_bits
            .as_ref()
            .is_some_and(|v| v.len() == value.atom_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn historical_morgan_recipes_retain_original_source_identity() {
        let value: serde_json::Value = serde_json::from_slice(
            &crate::read(&crate::directory().join("testdata/profiles/morgan_historical_0_3.json"))
                .unwrap(),
        )
        .unwrap();
        let recipes: Vec<Recipe> = serde_json::from_value(value["recipes"].clone()).unwrap();
        assert_eq!(recipes, Recipe::ALL);
        assert_eq!(
            value["source_sha256"].as_str().unwrap(),
            crate::digest(
                &crate::read(&crate::root().join("dev/tools/chembl_parity/audit_core.py")).unwrap()
            )
        );
    }
    #[test]
    fn historical_morgan_source_eligibility_and_public_calls() {
        for (smiles, expected) in [
            ("", 10),
            ("C", 14),
            ("CC", 15),
            ("C.C", 14),
            ("C[C@H](N)C(=O)O", 15),
            ("F/C=C/F", 15),
            ("c1ccccc1", 15),
        ] {
            let mol = Molecule::from_smiles(smiles).unwrap();
            let recipes: Vec<_> = Recipe::ALL
                .into_iter()
                .filter(|p| p.eligible(mol.num_atoms(), mol.num_bonds()))
                .collect();
            assert_eq!(recipes.len(), expected, "{smiles}");
            for recipe in recipes {
                let value = run(&mol, recipe).unwrap();
                assert!(valid(recipe, &value), "{smiles}:{recipe:?}");
            }
        }
        assert!(!Recipe::ChiralBondInvariants.params().include_chirality);
        assert_eq!(Recipe::Default.params().radius, 2);
        assert_eq!(Recipe::Chiral.params().radius, 2);
        assert_eq!(Recipe::Radius3.params().radius, 3);
        assert_eq!(Recipe::NoBondTypes.params().fp_size, 1024);
        assert!(!Recipe::NoBondTypes.params().use_bond_types);
        assert_eq!(Recipe::TwoBitsPerFeature.params().bits_per_feature, 2);
        assert!(
            Recipe::OnlyNonzeroCustomInvariants
                .params()
                .only_nonzero_invariants
        );
    }
}
