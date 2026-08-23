//! Family-neutral RDKit fingerprint-generator pipeline.
//!
//! Fingerprint families own invariant and environment construction. This module
//! owns the single environment-to-vector projection used by every family.

use super::{
    AdditionalOutput, Fingerprint, FingerprintArguments, FingerprintError,
    FingerprintFuncArguments, MorganArguments, MorganAtomInvGenerator,
    MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganFingerprintGenerator,
    RdkitFingerprintMtRng, SparseBitFingerprint, SparseCountFingerprint,
    atom_pair::{
        AtomPairArguments, AtomPairAtomInvariantsGenerator, AtomPairFingerprintGenerator,
        atom_pair_generator,
    },
    duplicate_additional_output_bit, getMorganGenerator, setup_temp_additional_output,
};

fn set_sparse_count_value(
    result: &mut SparseCountFingerprint,
    bit_id: u64,
    value: i32,
) -> Result<(), FingerprintError> {
    // RDKit source file: DataStructs/SparseIntVect.h
    // RDKit source: SparseIntVect.h lines 90-99
    // RDKit✔️✔️: void setVal(IndexType idx, int val) {
    // RDKit✔️✔️:   if (!checkIndex(idx)) {
    if bit_id >= result.size() {
        // RDKit✔️✔️:     throw IndexErrorException(static_cast<int>(idx));
        return Err(FingerprintError::SparseIndexOutOfRange {
            index: bit_id,
            size: result.size(),
        });
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (val != 0) {
    // RDKit✔️✔️:     d_data[idx] = val;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     d_data.erase(idx);
    // RDKit✔️✔️:   }
    result.set_val(bit_id, value);
    // RDKit✔️✔️: }
    Ok(())
}
use crate::{Molecule, notation::smiles_write::assign_stereochemistry_on_working_copy};
use serde_json::Value;

/// One source fingerprint environment after family-specific enumeration.
pub(crate) trait FingerprintEnvironment {
    fn bit_id(
        &self,
        arguments: &FingerprintArguments,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        hash_results: bool,
        fp_size: u64,
    ) -> Result<u64, FingerprintError>;

    fn update_additional_output(&self, output: &mut AdditionalOutput, bit_id: u64);
}

/// Family-specific work consumed by the common RDKit generator pipeline.
pub(crate) trait FingerprintFamily {
    type Environment: FingerprintEnvironment;

    fn common_arguments(&self) -> &FingerprintArguments;
    fn result_size(&self) -> u64;
    fn arguments_info_string(&self) -> String;
    fn environment_info_string(&self) -> String;
    fn atom_invariants_info_string(&self) -> Option<String>;
    fn bond_invariants_info_string(&self) -> Option<String>;
    fn arguments_json(&self) -> String;
    fn environment_json(&self) -> String;
    fn atom_invariants_json(&self) -> Option<String>;
    fn bond_invariants_json(&self) -> Option<String>;

    fn atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError>;
    fn bond_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError>;

    #[allow(clippy::too_many_arguments)]
    fn environments(
        &self,
        molecule: &Molecule,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        conf_id: i32,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        hash_results: bool,
    ) -> Result<Vec<Self::Environment>, FingerprintError>;
}

/// The sole common fingerprint result projector.
pub(crate) struct FingerprintGenerator<'a, F: FingerprintFamily> {
    family: &'a F,
}

#[derive(Debug, Clone)]
pub(crate) enum RestoredFingerprintGenerator {
    Morgan(MorganFingerprintGenerator),
    AtomPair(AtomPairFingerprintGenerator),
}

impl RestoredFingerprintGenerator {
    pub(crate) fn info_string(&self) -> String {
        match self {
            Self::Morgan(generator) => FingerprintGenerator::new(generator).info_string(),
            Self::AtomPair(generator) => generator.info_string(),
        }
    }

    pub(crate) fn to_json(&self) -> String {
        match self {
            Self::Morgan(generator) => FingerprintGenerator::new(generator).to_json(),
            Self::AtomPair(generator) => generator.to_json(),
        }
    }

    pub(crate) fn sparse_count_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        match self {
            Self::Morgan(generator) => generator.getSparseCountFingerprint(molecule, arguments),
            Self::AtomPair(generator) => generator.sparse_count_fingerprint(molecule, arguments),
        }
    }

    pub(crate) fn sparse_bit_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseBitFingerprint, FingerprintError> {
        match self {
            Self::Morgan(generator) => generator.getSparseFingerprint(molecule, arguments),
            Self::AtomPair(generator) => generator.sparse_bit_fingerprint(molecule, arguments),
        }
    }

    pub(crate) fn folded_count_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        match self {
            Self::Morgan(generator) => generator.getCountFingerprint(molecule, arguments),
            Self::AtomPair(generator) => generator.count_fingerprint(molecule, arguments),
        }
    }

    pub(crate) fn explicit_bit_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<Fingerprint, FingerprintError> {
        match self {
            Self::Morgan(generator) => generator.getFingerprint(molecule, arguments),
            Self::AtomPair(generator) => generator.fingerprint(molecule, arguments),
        }
    }
}

fn required_generator_child<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
) -> Result<&'a serde_json::Map<String, Value>, FingerprintError> {
    object.get(field).and_then(Value::as_object).ok_or_else(|| {
        FingerprintError::InvalidArgumentsJson(format!("{field} must be a JSON object"))
    })
}

fn required_generator_type<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
) -> Result<&'a str, FingerprintError> {
    object.get("type").and_then(Value::as_str).ok_or_else(|| {
        FingerprintError::InvalidArgumentsJson(format!("{field} type not specified in JSON"))
    })
}

pub(crate) fn generator_from_json(
    json: &str,
) -> Result<RestoredFingerprintGenerator, FingerprintError> {
    // RDKit source file: FingerprintGenerator.cpp
    // RDKit source: FingerprintGenerator.cpp lines 227-321
    // RDKit✔️✔️: std::unique_ptr<FingerprintGenerator<std::uint64_t>> generatorFromJSON(
    // RDKit✔️✔️:     const std::string &json) {
    // RDKit✔️✔️:   std::istringstream ss;
    // RDKit✔️✔️:   ss.str(json);
    // RDKit✔️✔️:   boost::property_tree::ptree pt;
    // RDKit✔️✔️:   boost::property_tree::read_json(ss, pt);
    let value: Value = serde_json::from_str(json)
        .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
    let object = value.as_object().ok_or_else(|| {
        FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
    })?;

    // RDKit✔️✔️:   auto fpArgsNode = pt.get_child_optional("fingerprintArguments");
    // RDKit✔️✔️:   if (fpArgsNode) {
    let arguments_node = required_generator_child(object, "fingerprintArguments")?;
    // RDKit✔️✔️:     auto typ = fpArgsNode->get_optional<std::string>("type");
    // RDKit✔️✔️:     if (!typ) {
    // RDKit✔️✔️:       throw ValueErrorException(
    // RDKit✔️✔️:           "FingerprintArguments type not specified in JSON");
    // RDKit✔️✔️:     }
    let arguments_type = required_generator_type(arguments_node, "FingerprintArguments")?;
    let arguments_json = Value::Object(arguments_node.clone()).to_string();

    // RDKit✔️✔️:   auto envGenNode = pt.get_child_optional("atomEnvironmentGenerator");
    // RDKit✔️✔️:   if (envGenNode) {
    let environment_node = required_generator_child(object, "atomEnvironmentGenerator")?;
    // RDKit✔️✔️:     auto typ = envGenNode->get_optional<std::string>("type");
    let environment_type = required_generator_type(environment_node, "AtomEnvironmentGenerator")?;

    // RDKit✔️✔️:   auto atomInvGenNode = pt.get_child_optional("atomInvariantsGenerator");
    let atom_invariants_node = object
        .get("atomInvariantsGenerator")
        .map(|value| {
            value.as_object().ok_or_else(|| {
                FingerprintError::InvalidArgumentsJson(
                    "atomInvariantsGenerator must be a JSON object".to_string(),
                )
            })
        })
        .transpose()?;

    match arguments_type {
        // RDKit✔️✔️:     if (*typ == "MorganArguments") {
        // RDKit✔️✔️:       fpArgs.reset(new MorganFingerprint::MorganArguments());
        "MorganArguments" => {
            if environment_type != "MorganEnvGenerator" {
                return Err(FingerprintError::InvalidArgumentsJson(format!(
                    "AtomEnvironmentGenerator type {environment_type} does not match MorganArguments"
                )));
            }
            let mut arguments = MorganArguments::default();
            arguments.fromJSON(&arguments_json)?;
            let atom_invariants_generator = match atom_invariants_node {
                None => None,
                Some(node) => match required_generator_type(node, "AtomInvariantsGenerator")? {
                    "MorganAtomInvGenerator" => {
                        let mut generator = MorganAtomInvGenerator::new(true);
                        generator.fromJSON(&Value::Object(node.clone()).to_string())?;
                        Some(MorganAtomInvariantsGenerator::Connectivity {
                            include_ring_membership: generator.include_ring_membership,
                        })
                    }
                    "MorganFeatureAtomInvGenerator" => {
                        if node.contains_key("patternSMARTS") {
                            return Err(FingerprintError::InvalidArgumentsJson(
                                "custom Morgan feature patterns are not representable by the current generator enum"
                                    .to_string(),
                            ));
                        }
                        Some(MorganAtomInvariantsGenerator::Feature)
                    }
                    unknown => {
                        return Err(FingerprintError::InvalidArgumentsJson(format!(
                            "Unknown AtomInvariantsGenerator type: {unknown}"
                        )));
                    }
                },
            };
            let bond_invariants_generator = object
                .get("bondInvariantsGenerator")
                .map(|value| {
                    let node = value.as_object().ok_or_else(|| {
                        FingerprintError::InvalidArgumentsJson(
                            "bondInvariantsGenerator must be a JSON object".to_string(),
                        )
                    })?;
                    let kind = required_generator_type(node, "BondInvariantsGenerator")?;
                    if kind != "MorganBondInvGenerator" {
                        return Err(FingerprintError::InvalidArgumentsJson(format!(
                            "Unknown BondInvariantsGenerator type: {kind}"
                        )));
                    }
                    let mut generator = MorganBondInvariantsGenerator::new(true, false);
                    generator.fromJSON(&Value::Object(node.clone()).to_string())?;
                    Ok(generator)
                })
                .transpose()?;
            Ok(RestoredFingerprintGenerator::Morgan(getMorganGenerator(
                &arguments,
                atom_invariants_generator,
                bond_invariants_generator,
                true,
                true,
            )))
        }
        // RDKit✔️✔️:     } else if (*typ == "AtomPairArguments") {
        // RDKit✔️✔️:       fpArgs.reset(new AtomPair::AtomPairArguments());
        "AtomPairArguments" => {
            if environment_type != "AtomPairEnvGenerator" {
                return Err(FingerprintError::InvalidArgumentsJson(format!(
                    "AtomEnvironmentGenerator type {environment_type} does not match AtomPairArguments"
                )));
            }
            if object.contains_key("bondInvariantsGenerator") {
                return Err(FingerprintError::InvalidArgumentsJson(
                    "AtomPair generator cannot contain a bondInvariantsGenerator".to_string(),
                ));
            }
            let mut arguments = AtomPairArguments::default();
            arguments.from_json(&arguments_json)?;
            let atom_invariants_generator = match atom_invariants_node {
                None => None,
                Some(node) => {
                    let kind = required_generator_type(node, "AtomInvariantsGenerator")?;
                    if kind != "AtomPairAtomInvGenerator" {
                        return Err(FingerprintError::InvalidArgumentsJson(format!(
                            "Unknown AtomInvariantsGenerator type: {kind}"
                        )));
                    }
                    let mut generator = AtomPairAtomInvariantsGenerator::default();
                    generator.from_json(&Value::Object(node.clone()).to_string())?;
                    Some(generator)
                }
            };
            Ok(RestoredFingerprintGenerator::AtomPair(atom_pair_generator(
                &arguments,
                atom_invariants_generator,
                true,
            )))
        }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       throw ValueErrorException("Unknown FingerprintArguments type: " + *typ);
        unknown => Err(FingerprintError::InvalidArgumentsJson(format!(
            "Unknown FingerprintArguments type: {unknown}"
        ))),
    }
    // RDKit✔️✔️: }
}

impl<'a, F: FingerprintFamily> FingerprintGenerator<'a, F> {
    pub(crate) fn new(family: &'a F) -> Self {
        Self { family }
    }

    pub(crate) fn info_string(&self) -> String {
        // RDKit source file: FingerprintGenerator.cpp
        // RDKit source: FingerprintGenerator.cpp lines 162-173
        // RDKit✔️✔️: std::string FingerprintGenerator<OutputType>::infoString() const {
        // RDKit✔️✔️:   std::string separator = " --- ";
        let separator = " --- ";
        // RDKit✔️✔️:   return dp_fingerprintArguments->commonArgumentsString() + separator +
        // RDKit✔️✔️:          dp_fingerprintArguments->infoString() + separator +
        // RDKit✔️✔️:          dp_atomEnvironmentGenerator->infoString() + separator +
        // RDKit✔️✔️:          (dp_atomInvariantsGenerator
        // RDKit✔️✔️:               ? (dp_atomInvariantsGenerator->infoString() + separator)
        // RDKit✔️✔️:               : ("No atom invariants generator" + separator)) +
        // RDKit✔️✔️:          (dp_bondInvariantsGenerator
        // RDKit✔️✔️:               ? (dp_bondInvariantsGenerator->infoString())
        // RDKit✔️✔️:               : "No bond invariants generator");
        // RDKit✔️✔️: }
        format!(
            "{}{}{}{}{}{}{}{}{}",
            self.family.common_arguments().common_arguments_string(),
            separator,
            self.family.arguments_info_string(),
            separator,
            self.family.environment_info_string(),
            separator,
            self.family
                .atom_invariants_info_string()
                .unwrap_or_else(|| "No atom invariants generator".to_string()),
            separator,
            self.family
                .bond_invariants_info_string()
                .unwrap_or_else(|| "No bond invariants generator".to_string())
        )
    }

    pub(crate) fn to_json(&self) -> String {
        // RDKit source: FingerprintGenerator.cpp lines 181-201
        // RDKit✔️✔️: void FingerprintGenerator<OutputType>::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("name", "FingerprintGenerator");
        let mut fields = vec![r#""name":"FingerprintGenerator""#.to_string()];
        // RDKit✔️✔️:   boost::property_tree::ptree argsNode;
        // RDKit✔️✔️:   dp_fingerprintArguments->toJSON(argsNode);
        // RDKit✔️✔️:   pt.add_child("fingerprintArguments", argsNode);
        fields.push(format!(
            r#""fingerprintArguments":{}"#,
            self.family.arguments_json()
        ));
        // RDKit✔️✔️:   boost::property_tree::ptree envGenNode;
        // RDKit✔️✔️:   dp_atomEnvironmentGenerator->toJSON(envGenNode);
        // RDKit✔️✔️:   pt.add_child("atomEnvironmentGenerator", envGenNode);
        fields.push(format!(
            r#""atomEnvironmentGenerator":{}"#,
            self.family.environment_json()
        ));
        // RDKit✔️✔️:   if (dp_atomInvariantsGenerator) {
        // RDKit✔️✔️:     boost::property_tree::ptree atomInvGenNode;
        // RDKit✔️✔️:     dp_atomInvariantsGenerator->toJSON(atomInvGenNode);
        // RDKit✔️✔️:     pt.add_child("atomInvariantsGenerator", atomInvGenNode);
        // RDKit✔️✔️:   }
        if let Some(json) = self.family.atom_invariants_json() {
            fields.push(format!(r#""atomInvariantsGenerator":{json}"#));
        }
        // RDKit✔️✔️:   if (dp_bondInvariantsGenerator) {
        // RDKit✔️✔️:     boost::property_tree::ptree bondInvGenNode;
        // RDKit✔️✔️:     dp_bondInvariantsGenerator->toJSON(bondInvGenNode);
        // RDKit✔️✔️:     pt.add_child("bondInvariantsGenerator", bondInvGenNode);
        // RDKit✔️✔️:   }
        if let Some(json) = self.family.bond_invariants_json() {
            fields.push(format!(r#""bondInvariantsGenerator":{json}"#));
        }
        // RDKit✔️✔️: }
        format!("{{{}}}", fields.join(","))
    }

    pub(crate) fn from_json(&self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 203-209
        // RDKit✔️✔️: void FingerprintGenerator<OutputType>::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &) {}
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        Ok(())
    }

    pub(crate) fn get_sparse_count_fingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<OutputType>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getSparseCountFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   return getFingerprintHelper(mol, args);
        // RDKit✔️✔️: }
        self.get_fingerprint_helper(molecule, args, 0)
    }

    pub(crate) fn get_sparse_fingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseBitFingerprint, FingerprintError> {
        // RDKit✔️✔️: // todo getSparseFingerprint does not completely produce the same output as
        // RDKit✔️✔️: // getSparseCountFingerprint. Count simulation and potential 64 bit outputs
        // RDKit✔️✔️: // makes size limiting necessary for getSparseFingerprint. This can be
        // RDKit✔️✔️: // changed if there is another way to avoid the size limitation of
        // RDKit✔️✔️: // SparseBitVect
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseBitVect>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getSparseFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   // make sure the result will fit into SparseBitVect
        // RDKit✔️✔️:   std::uint32_t resultSize =
        // RDKit✔️✔️:       std::min((std::uint64_t)std::numeric_limits<std::uint32_t>::max(),
        // RDKit✔️✔️:                (std::uint64_t)dp_atomEnvironmentGenerator->getResultSize());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::uint32_t effectiveSize = resultSize;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:     // effective size needs to be smaller than result size to compansate for
        // RDKit✔️✔️:     // count simulation
        // RDKit✔️✔️:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   AdditionalOutput countSimulationOutput;
        // RDKit✔️✔️:   AdditionalOutput *origAO = nullptr;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
        // RDKit✔️✔️:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
        // RDKit✔️✔️:     origAO = args.additionalOutput;
        // RDKit✔️✔️:     args.additionalOutput = &countSimulationOutput;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   auto result = std::make_unique<SparseBitVect>(resultSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        // RDKit✔️✔️:     if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:       for (unsigned int i = 0;
        // RDKit✔️✔️:            i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // RDKit✔️✔️:         // for every bound in the d_countBounds in dp_fingerprintArguments,
        // RDKit✔️✔️:         // set a bit if the occurrence count is equal or higher than the bound
        // RDKit✔️✔️:         // for that bit
        // RDKit✔️✔️:         const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        // RDKit✔️✔️:         if (val.second >= static_cast<int>(bounds_count[i])) {
        // RDKit✔️✔️:           OutputType nBitId = val.first * bounds_count.size() + i;
        // RDKit✔️✔️:           result->setBit(nBitId);
        // RDKit✔️✔️:           if (args.additionalOutput) {
        // RDKit✔️✔️:             duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
        // RDKit✔️✔️:                                          static_cast<OutputType>(val.first),
        // RDKit✔️✔️:                                          nBitId);
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       result->setBit(val.first);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (origAO) {
        // RDKit✔️✔️:     if (origAO->atomCounts) {
        // RDKit✔️✔️:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     args.additionalOutput = origAO;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        let result_size = self.family.result_size().min(u64::from(u32::MAX)) as u32;
        let fp_args = self.family.common_arguments();
        let effective_size = if fp_args.df_count_simulation {
            result_size / fp_args.d_count_bounds.len() as u32
        } else {
            result_size
        };

        let mut original_additional_output = None;
        if fp_args.df_count_simulation && args.additional_output.is_some() {
            let mut count_simulation_output = AdditionalOutput::default();
            setup_temp_additional_output(args, &mut count_simulation_output, molecule.num_atoms());
            original_additional_output = args.additional_output.take();
            args.additional_output = Some(count_simulation_output);
        }

        let temp_result = self.get_fingerprint_helper(molecule, args, u64::from(effective_size))?;
        let mut result = SparseBitFingerprint::new(u64::from(result_size));

        for (&bit_id, &count) in temp_result.nonzero_elements() {
            if fp_args.df_count_simulation {
                for (idx, &bound) in fp_args.d_count_bounds.iter().enumerate() {
                    if count >= bound as i32 {
                        let new_bit_id = bit_id * fp_args.d_count_bounds.len() as u64 + idx as u64;
                        result.set_bit(new_bit_id);
                        if let (Some(temp_output), Some(orig_output)) = (
                            args.additional_output.as_ref(),
                            original_additional_output.as_mut(),
                        ) {
                            duplicate_additional_output_bit(
                                temp_output,
                                orig_output,
                                bit_id,
                                new_bit_id,
                            )?;
                        }
                    }
                }
            } else {
                result.set_bit(bit_id);
            }
        }

        restore_original_additional_output(args, original_additional_output);
        Ok(result)
    }

    pub(crate) fn get_count_fingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<std::uint32_t>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getCountFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   auto tempResult =
        // RDKit✔️✔️:       getFingerprintHelper(mol, args, dp_fingerprintArguments->d_fpSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   auto result = std::make_unique<SparseIntVect<std::uint32_t>>(
        // RDKit✔️✔️:       dp_fingerprintArguments->d_fpSize);
        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        // RDKit✔️✔️:     result->setVal(val.first, val.second);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        let fp_size = self.family.common_arguments().d_fp_size;
        let temp_result = self.get_fingerprint_helper(molecule, args, u64::from(fp_size))?;
        let mut result = SparseCountFingerprint::new(u64::from(fp_size));
        for (&bit_id, &count) in temp_result.nonzero_elements() {
            set_sparse_count_value(&mut result, bit_id, count)?;
        }
        Ok(result)
    }

    pub(crate) fn get_fingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<Fingerprint, FingerprintError> {
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<ExplicitBitVect>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_countBounds.empty()) {
        // RDKit✔️✔️:       throw ValueErrorException("Count bounds are empty");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_countBounds.size() >= effectiveSize) {
        // RDKit✔️✔️:       throw ValueErrorException("Count bounds size is >= fingerprint size");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // effective size needs to be smaller than result size to compensate for
        // RDKit✔️✔️:     // count simulation
        // RDKit✔️✔️:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   AdditionalOutput countSimulationOutput;
        // RDKit✔️✔️:   AdditionalOutput *origAO = nullptr;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
        // RDKit✔️✔️:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
        // RDKit✔️✔️:     origAO = args.additionalOutput;
        // RDKit✔️✔️:     args.additionalOutput = &countSimulationOutput;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   auto result =
        // RDKit✔️✔️:       std::make_unique<ExplicitBitVect>(dp_fingerprintArguments->d_fpSize);
        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        // RDKit✔️✔️:     if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:       for (unsigned int i = 0;
        // RDKit✔️✔️:            i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // RDKit✔️✔️:         // for every bound in the d_countBounds in dp_fingerprintArguments,
        // RDKit✔️✔️:         // set a bit if the occurrence count is equal or higher than the bound
        // RDKit✔️✔️:         // for that bit
        // RDKit✔️✔️:         const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        // RDKit✔️✔️:         if (val.second >= static_cast<int>(bounds_count[i])) {
        // RDKit✔️✔️:           OutputType nBitId = val.first * bounds_count.size() + i;
        // RDKit✔️✔️:           result->setBit(nBitId);
        // RDKit✔️✔️:           if (args.additionalOutput) {
        // RDKit✔️✔️:             duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
        // RDKit✔️✔️:                                          static_cast<OutputType>(val.first),
        // RDKit✔️✔️:                                          nBitId);
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       result->setBit(val.first);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (origAO) {
        // RDKit✔️✔️:     if (origAO->atomCounts) {
        // RDKit✔️✔️:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     args.additionalOutput = origAO;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        let fp_args = self.family.common_arguments();
        let mut effective_size = fp_args.d_fp_size;
        if fp_args.df_count_simulation {
            if fp_args.d_count_bounds.is_empty() {
                return Err(FingerprintError::InvalidArguments {
                    reason: "Count bounds are empty",
                });
            }
            if fp_args.d_count_bounds.len() >= effective_size as usize {
                return Err(FingerprintError::InvalidArguments {
                    reason: "Count bounds size is >= fingerprint size",
                });
            }
            effective_size /= fp_args.d_count_bounds.len() as u32;
        }

        let mut original_additional_output = None;
        if fp_args.df_count_simulation && args.additional_output.is_some() {
            let mut count_simulation_output = AdditionalOutput::default();
            setup_temp_additional_output(args, &mut count_simulation_output, molecule.num_atoms());
            original_additional_output = args.additional_output.take();
            args.additional_output = Some(count_simulation_output);
        }

        let temp_result = self.get_fingerprint_helper(molecule, args, u64::from(effective_size))?;
        let mut on_bits = Vec::new();
        for (&bit_id, &count) in temp_result.nonzero_elements() {
            if fp_args.df_count_simulation {
                for (idx, &bound) in fp_args.d_count_bounds.iter().enumerate() {
                    if count >= bound as i32 {
                        let new_bit_id = bit_id * fp_args.d_count_bounds.len() as u64 + idx as u64;
                        on_bits.push(new_bit_id as usize);
                        if let (Some(temp_output), Some(orig_output)) = (
                            args.additional_output.as_ref(),
                            original_additional_output.as_mut(),
                        ) {
                            duplicate_additional_output_bit(
                                temp_output,
                                orig_output,
                                bit_id,
                                new_bit_id,
                            )?;
                        }
                    }
                }
            } else {
                on_bits.push(bit_id as usize);
            }
        }

        restore_original_additional_output(args, original_additional_output);
        Ok(Fingerprint::from_on_bits(
            fp_args.d_fp_size as usize,
            on_bits,
        ))
    }

    pub(crate) fn get_fingerprint_helper(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
        fp_size: u64,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<OutputType>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getFingerprintHelper(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args,
        // RDKit✔️✔️:     const std::uint64_t fpSize) const {
        // RDKit✔️✔️:   const ROMol *lmol = &mol;
        // RDKit✔️✔️:   std::unique_ptr<ROMol> tmol;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_includeChirality &&
        // RDKit✔️✔️:       !mol.hasProp(common_properties::_StereochemDone)) {
        // RDKit✔️✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
        // RDKit✔️✔️:     MolOps::assignStereochemistry(*tmol);
        // RDKit✔️✔️:     lmol = tmol.get();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (args.additionalOutput) {
        // RDKit✔️✔️:     reinitAdditionalOutput(*args.additionalOutput, mol.getNumAtoms());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   bool hashResults = false;
        // RDKit✔️✔️:   if (fpSize != 0) {
        // RDKit✔️✔️:     hashResults = true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
        // RDKit✔️✔️:   if (args.customAtomInvariants) {
        // RDKit✔️✔️:     atomInvariants.reset(
        // RDKit✔️✔️:         new std::vector<std::uint32_t>(*args.customAtomInvariants));
        // RDKit✔️✔️:   } else if (dp_atomInvariantsGenerator) {
        // RDKit✔️✔️:     atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
        // RDKit✔️✔️:   if (args.customBondInvariants) {
        // RDKit✔️✔️:     bondInvariants.reset(
        // RDKit✔️✔️:         new std::vector<std::uint32_t>(*args.customBondInvariants));
        // RDKit✔️✔️:   } else if (dp_bondInvariantsGenerator) {
        // RDKit✔️✔️:     bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // create all atom environments that will generate the bit-ids that will
        // RDKit✔️✔️:   // make up the fingerprint
        // RDKit✔️✔️:   auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(
        // RDKit✔️✔️:       *lmol, dp_fingerprintArguments, args.fromAtoms, args.ignoreAtoms,
        // RDKit✔️✔️:       args.confId, args.additionalOutput, atomInvariants.get(),
        // RDKit✔️✔️:       bondInvariants.get(), hashResults);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // allocate the result
        // RDKit✔️✔️:   auto res = std::make_unique<SparseIntVect<OutputType>>(
        // RDKit✔️✔️:       fpSize ? fpSize : dp_atomEnvironmentGenerator->getResultSize());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // define a mersenne twister with customized parameters.
        // RDKit✔️✔️:   // The standard parameters (used to create boost::mt19937)
        // RDKit✔️✔️:   // result in an RNG that's much too computationally intensive
        // RDKit✔️✔️:   // to seed.
        // RDKit✔️✔️:   // These are the parameters that have been used for the RDKit fingerprint.
        // RDKit✔️✔️:   typedef boost::random::mersenne_twister<std::uint32_t, 32, 4, 2, 31,
        // RDKit✔️✔️:                                           0x9908b0df, 11, 7, 0x9d2c5680, 15,
        // RDKit✔️✔️:                                           0xefc60000, 18, 3346425566U>
        // RDKit✔️✔️:       rng_type;
        // RDKit✔️✔️:   typedef boost::uniform_int<> distrib_type;
        // RDKit✔️✔️:   typedef boost::variate_generator<rng_type &, distrib_type> source_type;
        // RDKit✔️✔️:   std::unique_ptr<rng_type> generator;
        // RDKit✔️✔️:   //
        // RDKit✔️✔️:   // if we generate arbitrarily sized ints then mod them down to the
        // RDKit✔️✔️:   // appropriate size, we can guarantee that a fingerprint of
        // RDKit✔️✔️:   // size x has the same bits set as one of size 2x that's been folded
        // RDKit✔️✔️:   // in half.  This is a nice guarantee to have.
        // RDKit✔️✔️:   //
        // RDKit✔️✔️:   std::unique_ptr<distrib_type> dist;
        // RDKit✔️✔️:   std::unique_ptr<source_type> randomSource;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
        // RDKit✔️✔️:     // we will only create the RNG if we're going to need it
        // RDKit✔️✔️:     generator.reset(new rng_type(42u));
        // RDKit✔️✔️:     dist.reset(new distrib_type(0, INT_MAX));
        // RDKit✔️✔️:     randomSource.reset(new source_type(*generator, *dist));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // iterate over every atom environment and generate bit-ids that will make
        // RDKit✔️✔️:   // up the fingerprint
        // RDKit✔️✔️:   for (const auto env : atomEnvironments) {
        // RDKit✔️✔️:     OutputType seed = env->getBitId(dp_fingerprintArguments,
        // RDKit✔️✔️:                                     atomInvariants.get(), bondInvariants.get(),
        // RDKit✔️✔️:                                     args.additionalOutput, hashResults, fpSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     auto bitId = seed;
        // RDKit✔️✔️:     if (fpSize != 0) {
        // RDKit✔️✔️:       bitId %= fpSize;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     res->setVal(bitId, res->getVal(bitId) + 1);
        // RDKit✔️✔️:     if (args.additionalOutput) {
        // RDKit✔️✔️:       env->updateAdditionalOutput(args.additionalOutput, bitId);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     // do the additional bits if required:
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
        // RDKit✔️✔️:       generator->seed(static_cast<rng_type::result_type>(seed));
        // RDKit✔️✔️:
        // RDKit✔️✔️:       for (boost::uint32_t bitN = 1;
        // RDKit✔️✔️:            bitN < dp_fingerprintArguments->d_numBitsPerFeature; ++bitN) {
        // RDKit✔️✔️:         bitId = (*randomSource)();
        // RDKit✔️✔️:         if (fpSize != 0) {
        // RDKit✔️✔️:           bitId %= fpSize;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         res->setVal(bitId, res->getVal(bitId) + 1);
        // RDKit✔️✔️:         if (args.additionalOutput) {
        // RDKit✔️✔️:           env->updateAdditionalOutput(args.additionalOutput, bitId);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     delete env;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        let fp_args = self.family.common_arguments();

        let mut prepared_molecule;
        let environment_molecule =
            if fp_args.df_include_chirality && molecule.prop("_StereochemDone").is_none() {
                prepared_molecule = molecule.clone();
                assign_stereochemistry_on_working_copy(&mut prepared_molecule, true).map_err(
                    |error| FingerprintError::StereoPreparation {
                        reason: error.to_string(),
                    },
                )?;
                &prepared_molecule
            } else {
                molecule
            };

        if let Some(additional_output) = args.additional_output.as_mut() {
            additional_output.reset_for_atom_count(molecule.num_atoms());
        }
        let hash_results = fp_size != 0;
        let atom_invariants = if let Some(custom) = &args.custom_atom_invariants {
            custom.clone()
        } else {
            self.family.atom_invariants(molecule)?
        };
        let bond_invariants = if let Some(custom) = &args.custom_bond_invariants {
            custom.clone()
        } else {
            self.family.bond_invariants(molecule)?
        };
        let environments = self.family.environments(
            environment_molecule,
            args.from_atoms.as_deref(),
            args.ignore_atoms.as_deref(),
            args.conf_id,
            &atom_invariants,
            &bond_invariants,
            hash_results,
        )?;

        let result_size = if fp_size == 0 {
            self.family.result_size()
        } else {
            fp_size
        };
        let mut result = SparseCountFingerprint::new(result_size);
        let mut random_source =
            (fp_args.d_num_bits_per_feature > 1).then(|| RdkitFingerprintMtRng::new(42));

        for environment in environments {
            let seed = environment.bit_id(
                fp_args,
                &atom_invariants,
                &bond_invariants,
                hash_results,
                fp_size,
            )?;
            let bit_id = if hash_results { seed % fp_size } else { seed };
            let next_value = result.get_val(bit_id) + 1;
            set_sparse_count_value(&mut result, bit_id, next_value)?;
            if let Some(output) = args.additional_output.as_mut() {
                environment.update_additional_output(output, bit_id);
            }

            if let Some(random_source) = random_source.as_mut() {
                random_source.seed(seed as u32);
                for _ in 1..fp_args.d_num_bits_per_feature {
                    let random_bit_id = u64::from(random_source.uniform_int_0_to_i32_max());
                    let random_bit_id = if hash_results {
                        random_bit_id % fp_size
                    } else {
                        random_bit_id
                    };
                    let next_value = result.get_val(random_bit_id) + 1;
                    set_sparse_count_value(&mut result, random_bit_id, next_value)?;
                    if let Some(output) = args.additional_output.as_mut() {
                        environment.update_additional_output(output, random_bit_id);
                    }
                }
            }
        }

        Ok(result)
    }
}

fn restore_original_additional_output(
    args: &mut FingerprintFuncArguments,
    original: Option<AdditionalOutput>,
) {
    if let Some(mut original) = original {
        if original.atom_counts.is_some() {
            original.atom_counts = args
                .additional_output
                .as_ref()
                .and_then(|output| output.atom_counts.clone());
        }
        args.additional_output = Some(original);
    }
}
