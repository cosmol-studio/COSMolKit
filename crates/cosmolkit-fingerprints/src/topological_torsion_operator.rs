//! Reusable detached operator state. ROOT decision fp-generator-shared-state.
use super::*;
use std::sync::{Arc, Mutex, MutexGuard};
#[derive(Debug)]
struct State {
    common: FingerprintArguments,
    torsion_atom_count: u32,
    only_shortest_paths: bool,
    atom_invariants: Option<AtomPairAtomInvariantsGenerator>,
}
impl State {
    fn params(&self) -> TopologicalTorsionParams {
        TopologicalTorsionParams {
            torsion_atom_count: self.torsion_atom_count,
            only_shortest_paths: self.only_shortest_paths,
            include_chirality: self.common.include_chirality,
            count_simulation: self.common.count_simulation,
            fp_size: self.common.fp_size,
            bits_per_feature: self.common.bits_per_feature,
            count_bounds: self.common.count_bounds.clone(),
        }
    }
    fn derived(&self) -> TopologicalTorsionParams {
        // Computation borrows canonical common state; this view contains only
        // primitive derived values. Empty local bounds are never used. No
        // per-call constructor, vector clone or option validation is introduced.
        TopologicalTorsionParams {
            torsion_atom_count: self.torsion_atom_count,
            only_shortest_paths: self.only_shortest_paths,
            include_chirality: self.common.include_chirality,
            count_bounds: Vec::new(),
            count_simulation: self.common.count_simulation,
            fp_size: self.common.fp_size,
            bits_per_feature: self.common.bits_per_feature,
        }
    }
}
/// Clone aliases the source-defined reusable generator state.
#[derive(Debug, Clone)]
pub struct TopologicalTorsionGenerator {
    state: Arc<Mutex<State>>,
}
/// Bound settings view; immutable detached call parameters stay separate.
#[derive(Debug, Clone)]
pub struct TopologicalTorsionSettings {
    state: Arc<Mutex<State>>,
}
fn lock(state: &Mutex<State>) -> Result<MutexGuard<'_, State>, TopologicalTorsionError> {
    state
        .lock()
        .map_err(|_| TopologicalTorsionError::StatePoisoned)
}
impl TopologicalTorsionGenerator {
    pub fn new(
        params: &TopologicalTorsionParams,
        atom_invariants: Option<AtomPairAtomInvariantsGenerator>,
    ) -> Result<Self, TopologicalTorsionError> {
        // Source constructor validation occurs once. getResultSize is not
        // called by RDKit's factory: e.g. torsionAtomCount=8 constructs, and
        // only its unfolded result-size calculation executes an undefined shift.
        // RDKit❗✔️: FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
        // RDKit❗✔️:     const TopologicalTorsionArguments &args,
        // RDKit❗✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
        // RDKit❗✔️:     const bool ownsAtomInvGen) {
        // RDKit❗✔️:   auto *envGenerator = new TopologicalTorsionEnvGenerator<OutputType>();
        // RDKit❗✔️:
        // RDKit❗✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
        // RDKit❗✔️:   if (!atomInvariantsGenerator) {
        // RDKit❗✔️:     atomInvariantsGenerator =
        // RDKit❗✔️:         new AtomPair::AtomPairAtomInvGenerator(args.df_includeChirality, true);
        // RDKit❗✔️:     ownsAtomInvGenerator = true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return new FingerprintGenerator<OutputType>(
        // RDKit❗✔️:       envGenerator, new TopologicalTorsionArguments(args),
        // RDKit❗✔️:       atomInvariantsGenerator, nullptr, ownsAtomInvGenerator, false);
        // RDKit❗✔️: };
        // RDKit❗✔️: template <typename OutputType>
        // RDKit❗✔️: FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
        // RDKit❗✔️:     bool includeChirality, uint32_t torsionAtomCount,
        // RDKit❗✔️:     AtomInvariantsGenerator *atomInvariantsGenerator, bool countSimulation,
        // RDKit❗✔️:     std::uint32_t fpSize, std::vector<std::uint32_t> countBounds,
        // RDKit❗✔️:     bool ownsAtomInvGen) {
        // RDKit❗✔️:   TopologicalTorsionArguments arguments(includeChirality, torsionAtomCount,
        // RDKit❗✔️:                                         countSimulation, countBounds, fpSize);
        // RDKit❗✔️:   return getTopologicalTorsionGenerator<OutputType>(
        // RDKit❗✔️:       arguments, atomInvariantsGenerator, ownsAtomInvGen);
        // RDKit❗✔️: };
        // RDKit❗✔️:
        // RDKit❗✔️: // Topological torsion fingerprint does not support 32 bit output yet
        // RDKit❗✔️:
        let common = params.common()?;
        let inv = atom_invariants.unwrap_or(AtomPairAtomInvariantsGenerator {
            include_chirality: params.include_chirality,
            topological_torsion_correction: true,
        });
        Ok(Self {
            state: Arc::new(Mutex::new(State {
                common,
                torsion_atom_count: params.torsion_atom_count,
                only_shortest_paths: params.only_shortest_paths,
                atom_invariants: Some(inv),
            })),
        })
    }
    pub fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings {
            state: Arc::clone(&self.state),
        }
    }
    pub fn info_string(&self) -> Result<String, TopologicalTorsionError> {
        // RDKit❗✔️: std::string FingerprintGenerator<OutputType>::infoString() const {
        // RDKit❗✔️:   std::string separator = " --- ";
        // RDKit❗✔️:   return dp_fingerprintArguments->commonArgumentsString() + separator +
        // RDKit❗✔️:          dp_fingerprintArguments->infoString() + separator +
        // RDKit❗✔️:          dp_atomEnvironmentGenerator->infoString() + separator +
        // RDKit❗✔️:          (dp_atomInvariantsGenerator
        // RDKit❗✔️:               ? (dp_atomInvariantsGenerator->infoString() + separator)
        // RDKit❗✔️:               : ("No atom invariants generator" + separator)) +
        // RDKit❗✔️:          (dp_bondInvariantsGenerator
        // RDKit❗✔️:               ? (dp_bondInvariantsGenerator->infoString())
        // RDKit❗✔️:               : "No bond invariants generator");
        // RDKit❗✔️: }
        // RDKit❗✔️:
        let state = lock(&self.state)?;
        let p = state.derived();
        Ok(format!(
            "{} --- {} --- TopologicalTorsionEnvGenerator --- {} --- No bond invariants generator",
            common_arguments_string(
                state.common.count_simulation,
                state.common.fp_size,
                state.common.bits_per_feature,
                state.common.include_chirality
            ),
            p.info_string(),
            state
                .atom_invariants
                .map(|i| i.info_string())
                .unwrap_or_else(|| "No atom invariants generator".into())
        ))
    }
    pub fn to_json(&self) -> Result<String, TopologicalTorsionError> {
        // RDKit❗✔️: void FingerprintGenerator<OutputType>::toJSON(
        // RDKit❗✔️:     boost::property_tree::ptree &pt) const {
        // RDKit❗✔️:   pt.put("name", "FingerprintGenerator");
        // RDKit❗✔️:   boost::property_tree::ptree argsNode;
        // RDKit❗✔️:   dp_fingerprintArguments->toJSON(argsNode);
        // RDKit❗✔️:   pt.add_child("fingerprintArguments", argsNode);
        // RDKit❗✔️:   boost::property_tree::ptree envGenNode;
        // RDKit❗✔️:   dp_atomEnvironmentGenerator->toJSON(envGenNode);
        // RDKit❗✔️:   pt.add_child("atomEnvironmentGenerator", envGenNode);
        // RDKit❗✔️:   if (dp_atomInvariantsGenerator) {
        // RDKit❗✔️:     boost::property_tree::ptree atomInvGenNode;
        // RDKit❗✔️:     dp_atomInvariantsGenerator->toJSON(atomInvGenNode);
        // RDKit❗✔️:     pt.add_child("atomInvariantsGenerator", atomInvGenNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (dp_bondInvariantsGenerator) {
        // RDKit❗✔️:     boost::property_tree::ptree bondInvGenNode;
        // RDKit❗✔️:     dp_bondInvariantsGenerator->toJSON(bondInvGenNode);
        // RDKit❗✔️:     pt.add_child("bondInvariantsGenerator", bondInvGenNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        let state = lock(&self.state)?;
        let mut result = format!(
            "{{\"name\":\"FingerprintGenerator\",\"fingerprintArguments\":{},\"atomEnvironmentGenerator\":{{\"type\":\"TopologicalTorsionEnvGenerator\"}}",
            state.params().to_json()
        );
        if let Some(inv) = state.atom_invariants {
            result.push_str(",\"atomInvariantsGenerator\":");
            result.push_str(&inv.to_json());
        }
        result.push('}');
        Ok(result)
    }
    pub fn from_json(json: &str) -> Result<Self, TopologicalTorsionError> {
        // RDKit❗✔️: std::unique_ptr<FingerprintGenerator<std::uint64_t>> generatorFromJSON(
        // RDKit❗✔️:     const std::string &json) {
        // RDKit❗✔️:   std::istringstream ss;
        // RDKit❗✔️:   ss.str(json);
        // RDKit❗✔️:   boost::property_tree::ptree pt;
        // RDKit❗✔️:   boost::property_tree::read_json(ss, pt);
        // RDKit❗✔️:
        // RDKit❗✔️:   std::unique_ptr<AtomEnvironmentGenerator<std::uint64_t>> envGen;
        // RDKit❗✔️:   std::unique_ptr<FingerprintArguments> fpArgs;
        // RDKit❗✔️:   std::unique_ptr<AtomInvariantsGenerator> atomInvGen;
        // RDKit❗✔️:   std::unique_ptr<BondInvariantsGenerator> bondInvGen;
        // RDKit❗✔️:
        // RDKit❗✔️:   auto fpArgsNode = pt.get_child_optional("fingerprintArguments");
        // RDKit❗✔️:   if (fpArgsNode) {
        // RDKit❗✔️:     auto typ = fpArgsNode->get_optional<std::string>("type");
        // RDKit❗✔️:     if (!typ) {
        // RDKit❗✔️:       throw ValueErrorException(
        // RDKit❗✔️:           "FingerprintArguments type not specified in JSON");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (*typ == "MorganArguments") {
        // RDKit❗✔️:       fpArgs.reset(new MorganFingerprint::MorganArguments());
        // RDKit❗✔️:     } else if (*typ == "RDKitFPArguments") {
        // RDKit❗✔️:       fpArgs.reset(new RDKitFP::RDKitFPArguments());
        // RDKit❗✔️:     } else if (*typ == "AtomPairArguments") {
        // RDKit❗✔️:       fpArgs.reset(new AtomPair::AtomPairArguments());
        // RDKit❗✔️:     } else if (*typ == "TopologicalTorsionArguments") {
        // RDKit❗✔️:       fpArgs.reset(new TopologicalTorsion::TopologicalTorsionArguments());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       throw ValueErrorException("Unknown FingerprintArguments type: " + *typ);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     fpArgs->fromJSON(*fpArgsNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto envGenNode = pt.get_child_optional("atomEnvironmentGenerator");
        // RDKit❗✔️:   if (envGenNode) {
        // RDKit❗✔️:     auto typ = envGenNode->get_optional<std::string>("type");
        // RDKit❗✔️:     if (!typ) {
        // RDKit❗✔️:       throw ValueErrorException(
        // RDKit❗✔️:           "AtomEnvironmentGenerator type not specified in JSON");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (*typ == "MorganEnvGenerator") {
        // RDKit❗✔️:       envGen.reset(new MorganFingerprint::MorganEnvGenerator<std::uint64_t>());
        // RDKit❗✔️:     } else if (*typ == "RDKitFPEnvGenerator") {
        // RDKit❗✔️:       envGen.reset(new RDKitFP::RDKitFPEnvGenerator<std::uint64_t>());
        // RDKit❗✔️:     } else if (*typ == "AtomPairEnvGenerator") {
        // RDKit❗✔️:       envGen.reset(new AtomPair::AtomPairEnvGenerator<std::uint64_t>());
        // RDKit❗✔️:     } else if (*typ == "TopologicalTorsionEnvGenerator") {
        // RDKit❗✔️:       envGen.reset(new TopologicalTorsion::TopologicalTorsionEnvGenerator<
        // RDKit❗✔️:                    std::uint64_t>());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       throw ValueErrorException("Unknown AtomEnvGenerator type: " + *typ);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     envGen->fromJSON(*envGenNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto atomInvGenNode = pt.get_child_optional("atomInvariantsGenerator");
        // RDKit❗✔️:   if (atomInvGenNode) {
        // RDKit❗✔️:     auto typ = atomInvGenNode->get_optional<std::string>("type");
        // RDKit❗✔️:     if (!typ) {
        // RDKit❗✔️:       throw ValueErrorException(
        // RDKit❗✔️:           "AtomInvariantsGenerator type not specified in JSON");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (*typ == "MorganAtomInvGenerator") {
        // RDKit❗✔️:       atomInvGen.reset(new MorganFingerprint::MorganAtomInvGenerator());
        // RDKit❗✔️:     } else if (*typ == "MorganFeatureAtomInvGenerator") {
        // RDKit❗✔️:       atomInvGen.reset(new MorganFingerprint::MorganFeatureAtomInvGenerator());
        // RDKit❗✔️:     } else if (*typ == "RDKitFPAtomInvGenerator") {
        // RDKit❗✔️:       atomInvGen.reset(new RDKitFP::RDKitFPAtomInvGenerator());
        // RDKit❗✔️:     } else if (*typ == "AtomPairAtomInvGenerator") {
        // RDKit❗✔️:       atomInvGen.reset(new AtomPair::AtomPairAtomInvGenerator());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       throw ValueErrorException("Unknown AtomInvariantsGenerator type: " +
        // RDKit❗✔️:                                 *typ);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     atomInvGen->fromJSON(*atomInvGenNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto bondInvGenNode = pt.get_child_optional("bondInvariantsGenerator");
        // RDKit❗✔️:   if (bondInvGenNode) {
        // RDKit❗✔️:     auto typ = bondInvGenNode->get_optional<std::string>("type");
        // RDKit❗✔️:     if (!typ) {
        // RDKit❗✔️:       throw ValueErrorException(
        // RDKit❗✔️:           "BondInvariantsGenerator type not specified in JSON");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (*typ == "MorganBondInvGenerator") {
        // RDKit❗✔️:       bondInvGen.reset(new MorganFingerprint::MorganBondInvGenerator());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       throw ValueErrorException("Unknown BondInvariantsGenerator type: " +
        // RDKit❗✔️:                                 *typ);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     bondInvGen->fromJSON(*bondInvGenNode);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return std::make_unique<FingerprintGenerator<std::uint64_t>>(
        // RDKit❗✔️:       envGen.release(), fpArgs.release(),
        // RDKit❗✔️:       atomInvGen ? atomInvGen.release() : nullptr,
        // RDKit❗✔️:       bondInvGen ? bondInvGen.release() : nullptr, true, true);
        // RDKit❗✔️: }
        // RDKit❗✔️:
        let value = parse_object(json)?;
        let child = |field: &str, category: &str| -> Result<&Value, TopologicalTorsionError> {
            value.get(field).filter(|v| v.is_object()).ok_or_else(|| {
                FingerprintJsonError::Invalid(format!("missing {category} node in JSON")).into()
            })
        };
        let args = child("fingerprintArguments", "FingerprintArguments")?;
        let env = child("atomEnvironmentGenerator", "AtomEnvironmentGenerator")?;
        let typ = |node: &Value, category: &str| -> Result<String, TopologicalTorsionError> {
            node.get("type")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .ok_or_else(|| {
                    FingerprintJsonError::Invalid(format!("{category} type not specified in JSON"))
                        .into()
                })
        };
        let arg_type = typ(args, "FingerprintArguments")?;
        if arg_type != "TopologicalTorsionArguments" {
            return Err(FingerprintJsonError::Invalid(
                "JSON does not describe a Topological Torsion fingerprint generator".into(),
            )
            .into());
        }
        let env_type = typ(env, "AtomEnvironmentGenerator")?;
        if env_type != "TopologicalTorsionEnvGenerator" {
            return Err(FingerprintJsonError::UnsupportedComponent {
                component: "atomEnvironmentGenerator",
                source_type: env_type,
            }
            .into());
        }
        let mut params = TopologicalTorsionParams::default();
        params.from_json(&args.to_string())?;
        let atom_invariants = if let Some(node) = value.get("atomInvariantsGenerator") {
            let inv_type = typ(node, "AtomInvariantsGenerator")?;
            if inv_type != "AtomPairAtomInvGenerator" {
                return Err(FingerprintJsonError::UnsupportedComponent {
                    component: "atomInvariantsGenerator",
                    source_type: inv_type,
                }
                .into());
            }
            let mut inv = AtomPairAtomInvariantsGenerator::default();
            inv.from_json(&node.to_string())?;
            Some(inv)
        } else {
            None
        };
        if let Some(node) = value.get("bondInvariantsGenerator") {
            return Err(FingerprintJsonError::UnsupportedComponent {
                component: "bondInvariantsGenerator",
                source_type: typ(node, "BondInvariantsGenerator")?,
            }
            .into());
        }
        // Restoration mutates source default arguments without construction
        // preconditions; empty bounds and zero live options remain observable.
        Ok(Self {
            state: Arc::new(Mutex::new(State {
                common: FingerprintArguments {
                    count_simulation: params.count_simulation,
                    include_chirality: params.include_chirality,
                    count_bounds: params.count_bounds,
                    fp_size: params.fp_size,
                    bits_per_feature: params.bits_per_feature,
                },
                torsion_atom_count: params.torsion_atom_count,
                only_shortest_paths: params.only_shortest_paths,
                atom_invariants,
            })),
        })
    }
    pub fn sparse_count(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        Self::sparse_count_state(&state, input, call, output)
    }
    fn sparse_count_state(
        state: &State,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
        configured_count_helper(
            input,
            &state.derived(),
            &state.common,
            call,
            0,
            output,
            TorsionCodeMode::Modern,
            Some(state.atom_invariants),
        )
    }
    pub fn sparse_bits(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        Self::sparse_bits_state(&state, input, call, output)
    }
    fn sparse_bits_state(
        state: &State,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionError> {
        let p = state.derived();
        project_sparse_fingerprint(
            &state.common,
            p.result_size()?,
            input.topology.atoms.len(),
            output,
            |size, out| {
                configured_count_helper(
                    input,
                    &p,
                    &state.common,
                    call,
                    size,
                    out,
                    TorsionCodeMode::Modern,
                    Some(state.atom_invariants),
                )
            },
        )
    }
    pub fn count(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        Self::count_state(&state, input, call, output)
    }
    fn count_state(
        state: &State,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionError> {
        let p = state.derived();
        project_count_fingerprint(
            &state.common,
            input.topology.atoms.len(),
            output,
            |size, out| {
                configured_count_helper(
                    input,
                    &p,
                    &state.common,
                    call,
                    size,
                    out,
                    TorsionCodeMode::Modern,
                    Some(state.atom_invariants),
                )
            },
        )
    }
    pub fn bits(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        Self::bits_state(&state, input, call, output)
    }
    fn bits_state(
        state: &State,
        input: &AtomPairPreparedInput<'_>,
        call: &TopologicalTorsionCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, TopologicalTorsionError> {
        let p = state.derived();
        project_fingerprint(
            &state.common,
            input.topology.atoms.len(),
            output,
            |size, out| {
                configured_count_helper(
                    input,
                    &p,
                    &state.common,
                    call,
                    size,
                    out,
                    TorsionCodeMode::Modern,
                    Some(state.atom_invariants),
                )
            },
        )
    }
}
macro_rules! bound_setting {
    ($get:ident,$set:ident,$type:ty,$($field:ident).+) => {
        pub fn $get(&self)->Result<$type,TopologicalTorsionError> {Ok(lock(&self.state)?.$($field).+)}
        pub fn $set(&self,value:$type)->Result<(),TopologicalTorsionError> {lock(&self.state)?.$($field).+ =value;Ok(())}
    };
}
impl TopologicalTorsionSettings {
    pub fn snapshot(&self) -> Result<TopologicalTorsionParams, TopologicalTorsionError> {
        Ok(lock(&self.state)?.params())
    }
    bound_setting!(
        torsion_atom_count,
        set_torsion_atom_count,
        u32,
        torsion_atom_count
    );
    bound_setting!(
        only_shortest_paths,
        set_only_shortest_paths,
        bool,
        only_shortest_paths
    );
    bound_setting!(
        include_chirality,
        set_include_chirality,
        bool,
        common.include_chirality
    );
    bound_setting!(
        count_simulation,
        set_count_simulation,
        bool,
        common.count_simulation
    );
    bound_setting!(fp_size, set_fp_size, u32, common.fp_size);
    bound_setting!(
        bits_per_feature,
        set_bits_per_feature,
        u32,
        common.bits_per_feature
    );
    pub fn count_bounds(&self) -> Result<Vec<u32>, TopologicalTorsionError> {
        Ok(lock(&self.state)?.common.count_bounds.clone())
    }
    pub fn set_count_bounds(&self, values: Vec<u32>) -> Result<(), TopologicalTorsionError> {
        lock(&self.state)?.common.count_bounds = values;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestMolecule;
    #[test]
    fn original_generator_all_four_scalar_vector_forms() {
        let m = TestMolecule::from_smiles("CCCCO").unwrap();
        let input = m.input();
        let call = TopologicalTorsionCall::default();
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        let sparse = g.sparse_count(&input, &call, None).unwrap();
        assert_eq!(sparse.nonzero_elements().values().sum::<i32>(), 2);
        assert_eq!(
            sparse
                .nonzero_elements()
                .keys()
                .copied()
                .collect::<Vec<_>>(),
            vec![4_437_590_048, 12_893_306_913]
        );
        let sparse_bit = g.sparse_bits(&input, &call, None).unwrap();
        assert_eq!(sparse_bit.n_bits(), u32::MAX);
        assert_eq!(sparse_bit.on_bits().len(), 2);
        let bit = g.bits(&input, &call, None).unwrap();
        assert_eq!(bit.n_bits(), 2048);
        assert_eq!(bit.on_bits().len(), 2);
        let folded = TopologicalTorsionGenerator::new(
            &TopologicalTorsionParams {
                fp_size: 1000,
                ..Default::default()
            },
            None,
        )
        .unwrap()
        .count(&input, &call, None)
        .unwrap();
        assert_eq!(folded.length(), 1000);
        assert_eq!(folded.nonzero_elements().values().sum::<i32>(), 2);
        assert_eq!(
            folded
                .nonzero_elements()
                .keys()
                .copied()
                .collect::<Vec<_>>(),
            vec![24, 288]
        );
    }
    #[test]
    fn original_mutable_options_alias_canonical_generator_configuration() {
        let m = TestMolecule::from_smiles("CCCCO.Cl").unwrap();
        let input = m.input();
        let call = TopologicalTorsionCall::default();
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        let alias = g.clone();
        let settings = g.settings();
        let four = g.sparse_count(&input, &call, None).unwrap();
        assert_eq!(four.nonzero_elements().values().sum::<i32>(), 2);
        assert_eq!(four.length(), 1_u64 << 36);
        settings.set_torsion_atom_count(3).unwrap();
        settings.set_only_shortest_paths(true).unwrap();
        settings.set_fp_size(1024).unwrap();
        let three = alias.sparse_count(&input, &call, None).unwrap();
        let folded = alias.count(&input, &call, None).unwrap();
        assert_eq!(three.nonzero_elements().values().sum::<i32>(), 3);
        assert_eq!(three.length(), 1_u64 << 27);
        assert_ne!(three, four);
        assert_eq!(folded.length(), 1024);
        assert_eq!(folded.nonzero_elements().values().sum::<i32>(), 3);
        settings.set_include_chirality(true).unwrap();
        assert_eq!(
            g.sparse_count(&input, &call, None).unwrap().length(),
            1_u64 << 33
        );
        // Source factory creates atom generator once; mutating arguments does
        // not mutate an independent atom invariants generator's chirality.
        assert!(
            !lock(&g.state)
                .unwrap()
                .atom_invariants
                .unwrap()
                .include_chirality
        );
    }
    #[test]
    fn original_live_invalid_options_keep_valid_unfolded_counts() {
        let m = TestMolecule::from_smiles("CCCC").unwrap();
        let input = m.input();
        let call = TopologicalTorsionCall::default();
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        let settings = g.settings();
        settings.set_count_bounds(Vec::new()).unwrap();
        assert!(g.sparse_count(&input, &call, None).is_ok());
        assert!(matches!(
            g.bits(&input, &call, None),
            Err(TopologicalTorsionError::Fingerprint(
                FingerprintError::InvalidArguments {
                    reason: "Count bounds are empty"
                }
            ))
        ));
        // Dense source throws; sparse source executes undefined division.
        // Preserve the difference, using existing explicit safety category.
        assert!(matches!(
            g.sparse_bits(&input, &call, None),
            Err(TopologicalTorsionError::Fingerprint(
                FingerprintError::UndefinedArithmetic {
                    site: "FingerprintGenerator::getSparseFingerprint effectiveSize / countBounds.size()"
                }
            ))
        ));
        settings.set_count_simulation(false).unwrap();
        settings.set_fp_size(0).unwrap();
        assert!(g.sparse_count(&input, &call, None).is_ok());
        assert!(g.sparse_bits(&input, &call, None).is_ok());
        assert!(g.count(&input, &call, None).is_err());
        assert!(g.bits(&input, &call, None).is_err());
        settings.set_torsion_atom_count(8).unwrap();
        assert!(matches!(
            g.sparse_count(&input, &call, None),
            Err(TopologicalTorsionError::Fingerprint(
                FingerprintError::InvalidArguments {
                    reason: "topological torsion result-size shift must be less than 64 bits"
                }
            ))
        ));
    }
    #[test]
    fn original_json_roundtrip_retains_source_configuration_and_output() {
        let m = TestMolecule::from_smiles("C1CC1").unwrap();
        let params = TopologicalTorsionParams {
            include_chirality: true,
            torsion_atom_count: 3,
            count_bounds: vec![1, 3, 5],
            fp_size: 1536,
            only_shortest_paths: true,
            bits_per_feature: 2,
            ..Default::default()
        };
        let g = TopologicalTorsionGenerator::new(&params, None).unwrap();
        let json = g.to_json().unwrap();
        let restored = TopologicalTorsionGenerator::from_json(&json).unwrap();
        assert_eq!(restored.settings().snapshot().unwrap(), params);
        assert_eq!(
            lock(&restored.state).unwrap().atom_invariants,
            Some(AtomPairAtomInvariantsGenerator {
                include_chirality: true,
                topological_torsion_correction: true
            })
        );
        assert_eq!(
            restored
                .bits(&m.input(), &TopologicalTorsionCall::default(), None)
                .unwrap(),
            g.bits(&m.input(), &TopologicalTorsionCall::default(), None)
                .unwrap()
        );
    }
    #[test]
    fn original_supplied_atom_invariant_generator_is_retained() {
        let supplied = AtomPairAtomInvariantsGenerator {
            include_chirality: false,
            topological_torsion_correction: false,
        };
        let params = TopologicalTorsionParams::default();
        let g = TopologicalTorsionGenerator::new(&params, Some(supplied)).unwrap();
        assert_eq!(lock(&g.state).unwrap().atom_invariants, Some(supplied));
        let default = TopologicalTorsionGenerator::new(&params, None).unwrap();
        let m = TestMolecule::from_smiles("CCCCO").unwrap();
        assert_ne!(
            g.sparse_count(&m.input(), &TopologicalTorsionCall::default(), None)
                .unwrap(),
            default
                .sparse_count(&m.input(), &TopologicalTorsionCall::default(), None)
                .unwrap()
        );
    }
    #[test]
    fn source_native_metadata_and_absent_generator_restoration_have_no_fallback() {
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        assert_eq!(
            g.info_string().unwrap(),
            "Common arguments : countSimulation=1 fpSize=2048 bitsPerFeature=1 includeChirality=0 --- TopologicalTorsionArguments torsionAtomCount=4 onlyShortestPaths=0 --- TopologicalTorsionEnvGenerator --- AtomPairInvariantGenerator topologicalTorsionCorrection=1 --- No bond invariants generator"
        );
        let mut value: Value = serde_json::from_str(&g.to_json().unwrap()).unwrap();
        value
            .as_object_mut()
            .unwrap()
            .remove("atomInvariantsGenerator");
        value["fingerprintArguments"]["countBounds"] = Value::String(String::new());
        let restored = TopologicalTorsionGenerator::from_json(&value.to_string()).unwrap();
        assert!(
            restored
                .info_string()
                .unwrap()
                .contains("No atom invariants generator")
        );
        assert!(
            !restored
                .to_json()
                .unwrap()
                .contains("atomInvariantsGenerator")
        );
        assert!(restored.settings().count_bounds().unwrap().is_empty());
        let m = TestMolecule::from_smiles("CCCC").unwrap();
        let call = TopologicalTorsionCall {
            custom_atom_invariants: Some(&[17, 18, 19, 20]),
            ..Default::default()
        };
        assert_eq!(
            restored
                .sparse_count(&m.input(), &call, None)
                .unwrap()
                .nonzero_elements()
                .values()
                .sum::<i32>(),
            1
        );
    }
}

// Source strided bulk scheduler. Inputs are detached borrowed values and
// missing rows remain None; shared canonical generator state is read once.
fn bulk_with_state<T, F>(
    state: &State,
    inputs: &[Option<AtomPairPreparedInput<'_>>],
    workers: std::num::NonZeroUsize,
    action: F,
) -> Result<Vec<Option<T>>, TopologicalTorsionError>
where
    T: Send,
    F: Fn(
            &State,
            &AtomPairPreparedInput<'_>,
            &TopologicalTorsionCall<'_>,
        ) -> Result<T, TopologicalTorsionError>
        + Sync,
{
    // RDKit❗✔️: template <typename ReturnType, typename FuncType>
    // RDKit❗✔️: std::vector<std::unique_ptr<ReturnType>> mtgetFingerprints(
    // RDKit❗✔️:     FuncType func, const std::vector<const ROMol *> &mols, int numThreads) {
    // RDKit❗✔️:   std::vector<std::uint32_t> *fromAtoms = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *ignoreAtoms = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *customAtomInvariants = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *customBondInvariants = nullptr;
    // RDKit❗✔️:   int confId = -1;
    // RDKit❗✔️:   AdditionalOutput *additionalOutput = nullptr;
    // RDKit❗✔️:   FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
    // RDKit❗✔️:                                 additionalOutput, customAtomInvariants,
    // RDKit❗✔️:                                 customBondInvariants);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<std::unique_ptr<ReturnType>> result;
    // RDKit❗✔️:   auto numThreadsToUse = getNumThreadsToUse(numThreads);
    // RDKit❗✔️:   unsigned int nmols = mols.size();
    // RDKit❗✔️:   result.reserve(nmols);
    // RDKit❗✔️:   if (numThreadsToUse == 1) {
    // RDKit❗✔️:     for (auto i = 0u; i < nmols; ++i) {
    // RDKit❗✔️:       if (!mols[i]) {
    // RDKit❗✔️:         result.emplace_back(std::unique_ptr<ReturnType>());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         result.emplace_back(std::move(func(*mols[i], args)));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   else {
    // RDKit❗✔️:     std::vector<std::vector<std::unique_ptr<ReturnType>>> accum(
    // RDKit❗✔️:         numThreadsToUse);
    // RDKit❗✔️:     std::vector<std::thread> tg;
    // RDKit❗✔️:     for (auto ti = 0u; ti < numThreadsToUse; ++ti) {
    // RDKit❗✔️:       auto lfunc = [&](unsigned int tidx) {
    // RDKit❗✔️:         for (auto midx = tidx; midx < mols.size(); midx += numThreadsToUse) {
    // RDKit❗✔️:           if (!mols[midx]) {
    // RDKit❗✔️:             accum[tidx].emplace_back(std::unique_ptr<ReturnType>());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             accum[tidx].emplace_back(std::move(func(*mols[midx], args)));
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       };
    // RDKit❗✔️:       tg.emplace_back(std::thread(lfunc, ti));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto &thread : tg) {
    // RDKit❗✔️:       if (thread.joinable()) {
    // RDKit❗✔️:         thread.join();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto midx = 0u; midx < mols.size(); ++midx) {
    // RDKit❗✔️:       auto tidx = midx % numThreadsToUse;
    // RDKit❗✔️:       auto jidx = midx / numThreadsToUse;
    // RDKit❗✔️:       result.emplace_back(std::move(accum[tidx][jidx]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // Behavior: source ti/midx strides, every None slot, and output order
    // are retained. Caller already holds the sole state read lock. Each
    // worker uses the same configured generator/default call without
    // reconstructing/revalidating it. Errors stay structured; all launched
    // workers are joined before a spawn/calculation/panic error is returned.
    // Source C++ uncaught worker exceptions terminate its process; structured
    // worker errors are a separately reported canonical safety boundary.
    // Complexity: O(N) output movement and O(W) thread/iterator bookkeeping;
    // no Molecule/domain-state cloning or second algorithm path. Source-like
    // local contiguous vectors retain per-thread strided accumulation.
    let n = workers.get();
    let call = TopologicalTorsionCall::default();
    if n == 1 {
        return inputs
            .iter()
            .map(|input| {
                input
                    .as_ref()
                    .map(|input| action(state, input, &call))
                    .transpose()
            })
            .collect();
    }
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(n);
        let mut launch_error = None;
        for ti in 0..n {
            let action = &action;
            let call = &call;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                (ti..inputs.len())
                    .step_by(n)
                    .map(|midx| {
                        inputs[midx]
                            .as_ref()
                            .map(|input| action(state, input, call))
                            .transpose()
                    })
                    .collect::<Vec<_>>()
            }) {
                Ok(handle) => handles.push(handle),
                Err(error) => {
                    launch_error = Some(TopologicalTorsionError::ThreadSpawn(error));
                    break;
                }
            }
        }
        let mut accum = Vec::with_capacity(handles.len());
        let mut panic_error = false;
        for handle in handles {
            match handle.join() {
                Ok(rows) => accum.push(rows.into_iter()),
                Err(_) => panic_error = true,
            }
        }
        if let Some(error) = launch_error {
            return Err(error);
        }
        if panic_error {
            return Err(TopologicalTorsionError::WorkerPanic);
        }
        let mut result = Vec::with_capacity(inputs.len());
        for midx in 0..inputs.len() {
            let value = accum[midx % n]
                .next()
                .ok_or(TopologicalTorsionError::WorkerProtocol)?;
            result.push(value?);
        }
        Ok(result)
    })
}
impl TopologicalTorsionGenerator {
    pub fn sparse_counts(
        &self,
        inputs: &[Option<AtomPairPreparedInput<'_>>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseCountFingerprint>>, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        let workers = cosmolkit_core::rdkit_thread_count(num_threads)
            .map_err(TopologicalTorsionError::ThreadCount)?;
        let workers =
            std::num::NonZeroUsize::new(workers.get() as usize).expect("nonzero u32 worker count");
        bulk_with_state(&state, inputs, workers, |s, input, call| {
            Self::sparse_count_state(s, input, call, None)
        })
    }
    pub fn sparse_fingerprints(
        &self,
        inputs: &[Option<AtomPairPreparedInput<'_>>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseBitFingerprint>>, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        let workers = cosmolkit_core::rdkit_thread_count(num_threads)
            .map_err(TopologicalTorsionError::ThreadCount)?;
        let workers =
            std::num::NonZeroUsize::new(workers.get() as usize).expect("nonzero u32 worker count");
        bulk_with_state(&state, inputs, workers, |s, input, call| {
            Self::sparse_bits_state(s, input, call, None)
        })
    }
    pub fn counts(
        &self,
        inputs: &[Option<AtomPairPreparedInput<'_>>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseCountFingerprint32>>, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        let workers = cosmolkit_core::rdkit_thread_count(num_threads)
            .map_err(TopologicalTorsionError::ThreadCount)?;
        let workers =
            std::num::NonZeroUsize::new(workers.get() as usize).expect("nonzero u32 worker count");
        bulk_with_state(&state, inputs, workers, |s, input, call| {
            Self::count_state(s, input, call, None)
        })
    }
    pub fn fingerprints(
        &self,
        inputs: &[Option<AtomPairPreparedInput<'_>>],
        num_threads: i32,
    ) -> Result<Vec<Option<Fingerprint>>, TopologicalTorsionError> {
        let state = lock(&self.state)?;
        let workers = cosmolkit_core::rdkit_thread_count(num_threads)
            .map_err(TopologicalTorsionError::ThreadCount)?;
        let workers =
            std::num::NonZeroUsize::new(workers.get() as usize).expect("nonzero u32 worker count");
        bulk_with_state(&state, inputs, workers, |s, input, call| {
            Self::bits_state(s, input, call, None)
        })
    }
}

#[cfg(test)]
mod bulk_tests {
    use super::*;
    use crate::test_support::TestMolecule;
    #[test]
    fn source_strided_four_forms_keep_none_and_exact_input_order() {
        let a = TestMolecule::from_smiles("CCCCO").unwrap();
        let b = TestMolecule::from_smiles("CCCC").unwrap();
        let c = TestMolecule::from_smiles("C1CC1").unwrap();
        let inputs = [
            Some(a.input()),
            None,
            Some(b.input()),
            Some(c.input()),
            None,
        ];
        let generator =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        for workers in [1, 2, 7] {
            let rows = generator.sparse_counts(&inputs, workers).unwrap();
            assert_eq!(rows.len(), 5);
            assert!(rows[1].is_none() && rows[4].is_none());
            assert_eq!(
                rows[0]
                    .as_ref()
                    .unwrap()
                    .nonzero_elements()
                    .iter()
                    .map(|(&k, &v)| (k, v))
                    .collect::<Vec<_>>(),
                vec![(4437590048, 1), (12893306913, 1)]
            );
            assert_eq!(
                rows[2]
                    .as_ref()
                    .unwrap()
                    .nonzero_elements()
                    .iter()
                    .map(|(&k, &v)| (k, v))
                    .collect::<Vec<_>>(),
                vec![(4303372320, 1)]
            );
            assert_eq!(
                rows[3]
                    .as_ref()
                    .unwrap()
                    .nonzero_elements()
                    .iter()
                    .map(|(&k, &v)| (k, v))
                    .collect::<Vec<_>>(),
                vec![(4437590049, 1)]
            );
            for (m, result) in inputs
                .iter()
                .zip(generator.fingerprints(&inputs, workers).unwrap())
            {
                assert_eq!(
                    result,
                    m.as_ref().map(|i| generator
                        .bits(i, &TopologicalTorsionCall::default(), None)
                        .unwrap())
                );
            }
            for (m, result) in inputs
                .iter()
                .zip(generator.counts(&inputs, workers).unwrap())
            {
                assert_eq!(
                    result,
                    m.as_ref().map(|i| generator
                        .count(i, &TopologicalTorsionCall::default(), None)
                        .unwrap())
                );
            }
            for (m, result) in inputs
                .iter()
                .zip(generator.sparse_fingerprints(&inputs, workers).unwrap())
            {
                assert_eq!(
                    result,
                    m.as_ref().map(|i| generator
                        .sparse_bits(i, &TopologicalTorsionCall::default(), None)
                        .unwrap())
                );
            }
        }
    }
    #[test]
    fn source_empty_none_only_and_shared_live_options_bulk() {
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        assert!(g.fingerprints(&[], 2).unwrap().is_empty());
        assert_eq!(
            g.counts(&[None, None, None], 7).unwrap(),
            vec![None, None, None]
        );
        let a = TestMolecule::from_smiles("CCCCO").unwrap();
        let inputs = [Some(a.input()), None];
        g.settings().set_torsion_atom_count(3).unwrap();
        g.settings().set_fp_size(1000).unwrap();
        let rows = g.sparse_counts(&inputs, 2).unwrap();
        assert_eq!(rows[0].as_ref().unwrap().length(), 1 << 27);
        assert_eq!(
            rows[0]
                .as_ref()
                .unwrap()
                .nonzero_elements()
                .values()
                .sum::<i32>(),
            3
        );
        assert!(rows[1].is_none());
        assert_eq!(
            g.counts(&inputs, 2).unwrap()[0].as_ref().unwrap().length(),
            1000
        );
    }
    #[test]
    fn bulk_source_defined_dense_error_propagates_and_does_not_poison_state() {
        let g =
            TopologicalTorsionGenerator::new(&TopologicalTorsionParams::default(), None).unwrap();
        let a = TestMolecule::from_smiles("CCCC").unwrap();
        let inputs = [None, Some(a.input()), None];
        g.settings().set_count_bounds(vec![]).unwrap();
        assert!(matches!(
            g.fingerprints(&inputs, 2),
            Err(TopologicalTorsionError::Fingerprint(
                FingerprintError::InvalidArguments {
                    reason: "Count bounds are empty"
                }
            ))
        ));
        assert_eq!(
            g.sparse_counts(&inputs, 2).unwrap()[1]
                .as_ref()
                .unwrap()
                .nonzero_elements()
                .values()
                .sum::<i32>(),
            1
        );
        g.settings().set_count_bounds(vec![1, 2, 4, 8]).unwrap();
        assert!(g.fingerprints(&inputs, 2).is_ok());
        assert!(matches!(
            g.counts(&inputs, i32::MIN),
            Err(TopologicalTorsionError::ThreadCount(
                cosmolkit_core::ThreadCountError::UndefinedSignedNegation
            ))
        ));
    }
}
