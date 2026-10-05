//! Source-defined reusable Morgan generator state over the sole Morgan pipeline.
use super::*;
use crate::metadata::{SourceNode, bool_or, common_arguments_string, parse_object};
use crate::{AtomPairAtomInvariantsGenerator, AtomPairPreparedInput, FingerprintJsonError};
use cosmolkit_search::{
    QueryGraph, SmartsParseError, SmartsParseParams, SmartsWriteParams, parse_smarts,
    query_graph_to_smarts,
};
use serde_json::Value;
use std::sync::{Arc, Mutex, MutexGuard};

/// Independent immutable provider input; the factory captures it once.
#[derive(Debug, Clone)]
pub enum MorganAtomProvider {
    Connectivity { include_ring_membership: bool },
    Features { patterns: Option<Vec<QueryGraph>> },
    AtomPair(AtomPairAtomInvariantsGenerator),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MorganBondProvider {
    pub use_bond_types: bool,
    pub include_chirality: bool,
}
impl Default for MorganBondProvider {
    fn default() -> Self {
        Self {
            use_bond_types: true,
            include_chirality: false,
        }
    }
}
#[derive(Debug)]
enum AtomProvider {
    Connectivity,
    Features(MorganFeatureAtomInvGenerator),
    AtomPair(AtomPairAtomInvariantsGenerator),
}
#[derive(Debug)]
struct State {
    backend: MorganGenerator,
    atom_provider: Option<AtomProvider>,
    initial_ring_membership: bool,
    initial_use_bond_types: bool,
}
impl State {
    fn params(&self) -> MorganParams {
        let g = &self.backend;
        MorganParams {
            radius: g.radius,
            only_nonzero_invariants: g.only_nonzero_invariants,
            include_redundant_environments: g.include_redundant_environments,
            include_chirality: g.fingerprint_arguments.include_chirality,
            count_simulation: g.fingerprint_arguments.count_simulation,
            fp_size: g.fingerprint_arguments.fp_size,
            bits_per_feature: g.fingerprint_arguments.bits_per_feature,
            count_bounds: g.fingerprint_arguments.count_bounds.clone(),
            include_ring_membership: self.initial_ring_membership,
            use_bond_types: self.initial_use_bond_types,
        }
    }
    fn atoms(&self, input: &MorganPreparedInput<'_>) -> Result<Vec<u32>, MorganError> {
        // RDKit❗✔️:   } else if (dp_atomInvariantsGenerator) {
        // RDKit❗✔️:     atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
        // RDKit❗✔️:   }
        // Behavior: original molecule and prepared assignments, once-captured
        // independent provider; custom invariant slices bypass this callback.
        // Complexity: delegate sole invariant owners; no per-call query clone.
        match &self.atom_provider {
            None => Ok(Vec::new()),
            Some(AtomProvider::Connectivity) => self
                .backend
                .atom_invariants
                .as_ref()
                .expect("connectivity state has its configured owner")
                .get_atom_invariants(input.topology, input.valence, input.rings),
            Some(AtomProvider::Features(owner)) => {
                let target = SearchTarget::new(
                    input.topology,
                    input.coordinates,
                    &input.topology.stereo_groups,
                    Some(input.rings),
                    Some(input.valence),
                );
                let context =
                    build_prepared_query_match_context(input.topology, input.rings, input.valence)?;
                owner.get_atom_invariants(&target, &context)
            }
            Some(AtomProvider::AtomPair(owner)) => owner
                .atom_invariants(&AtomPairPreparedInput {
                    topology: input.topology,
                    properties: input.properties,
                    coordinates: input.coordinates,
                    valence: input.valence,
                    rings: input.rings,
                    use_legacy_stereo_perception: true,
                })
                .map_err(MorganError::from),
        }
    }
}
#[derive(Debug, Clone)]
pub struct MorganOperator {
    state: Arc<Mutex<State>>,
}
#[derive(Debug, Clone)]
pub struct MorganSettings {
    state: Arc<Mutex<State>>,
}
fn lock(state: &Mutex<State>) -> Result<MutexGuard<'_, State>, MorganError> {
    state.lock().map_err(|_| MorganError::StatePoisoned)
}

impl MorganOperator {
    pub fn new(
        params: &MorganParams,
        atoms: Option<MorganAtomProvider>,
        bonds: Option<MorganBondProvider>,
    ) -> Result<Self, MorganError> {
        // RDKit❗✔️:     const MorganArguments &args,
        // RDKit❗✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
        // RDKit❗✔️:     BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGen,
        // RDKit❗✔️:     bool ownsBondInvGen) {
        // RDKit❗✔️:   AtomEnvironmentGenerator<OutputType> *morganEnvGenerator =
        // RDKit❗✔️:       new MorganEnvGenerator<OutputType>();
        // RDKit❗✔️:
        // RDKit❗✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
        // RDKit❗✔️:   if (!atomInvariantsGenerator) {
        // RDKit❗✔️:     atomInvariantsGenerator = new MorganAtomInvGenerator();
        // RDKit❗✔️:     ownsAtomInvGenerator = true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   bool ownsBondInvGenerator = ownsBondInvGen;
        // RDKit❗✔️:   if (!bondInvariantsGenerator) {
        // RDKit❗✔️:     bondInvariantsGenerator = new MorganBondInvGenerator(
        // RDKit❗✔️:         args.df_useBondTypes, args.df_includeChirality);
        // RDKit❗✔️:     ownsBondInvGenerator = true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return new FingerprintGenerator<OutputType>(
        // RDKit❗✔️:       morganEnvGenerator, new MorganArguments(args), atomInvariantsGenerator,
        // RDKit❗✔️:       bondInvariantsGenerator, ownsAtomInvGenerator, ownsBondInvGenerator);
        // RDKit❗✔️: }
        // ROOT decision A: canonical include_ring_membership configures the
        // source explicit provider; bare Python wrapper's unused bool is separate.
        // Behavior: constructor common checks once, independent providers capture
        // initial options. Live option writes never reconstruct or refresh them.
        // Complexity: one common bounds copy and supplied pattern ownership copy;
        // subsequent calls borrow the persistent backend and provider owners.
        let mut backend = get_morgan_generator(params)?;
        let atom_provider = match atoms {
            None => AtomProvider::Connectivity,
            Some(MorganAtomProvider::Connectivity {
                include_ring_membership,
            }) => {
                backend.atom_invariants =
                    Some(MorganAtomInvGenerator::new(include_ring_membership));
                AtomProvider::Connectivity
            }
            Some(MorganAtomProvider::Features { patterns }) => {
                backend.atom_invariants = None;
                AtomProvider::Features(MorganFeatureAtomInvGenerator::from_owned_patterns(patterns))
            }
            Some(MorganAtomProvider::AtomPair(owner)) => {
                backend.atom_invariants = None;
                AtomProvider::AtomPair(owner)
            }
        };
        if let Some(p) = bonds {
            backend.bond_invariants = Some(MorganBondInvGenerator::new(
                p.use_bond_types,
                p.include_chirality,
            ));
        }
        Ok(Self {
            state: Arc::new(Mutex::new(State {
                backend,
                atom_provider: Some(atom_provider),
                initial_ring_membership: params.include_ring_membership,
                initial_use_bond_types: params.use_bond_types,
            })),
        })
    }
    pub fn settings(&self) -> MorganSettings {
        MorganSettings {
            state: Arc::clone(&self.state),
        }
    }
    pub fn info_string(&self) -> Result<String, MorganError> {
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
        // Behavior: full source sequence and separators, independent provider
        // metadata retained despite changes to live argument chirality.
        // Complexity: scalar/text formatting; no invariant calculation or clone.
        let state = lock(&self.state)?;
        let g = &state.backend;
        let common = &g.fingerprint_arguments;
        let atoms = match &state.atom_provider {
            None => "No atom invariants generator".into(),
            Some(AtomProvider::Connectivity) => format!(
                "MorganInvariantGenerator includeRingMembership={}",
                g.atom_invariants
                    .as_ref()
                    .expect("connectivity owner")
                    .include_ring_membership() as u8
            ),
            Some(AtomProvider::Features(_)) => "MorganFeatureInvariantGenerator".into(),
            Some(AtomProvider::AtomPair(owner)) => owner.info_string(),
        };
        let bonds = g
            .bond_invariants
            .as_ref()
            .map(|owner| {
                let (t, c) = owner.configuration();
                format!(
                    "MorganInvariantGenerator useBondTypes={} useChirality={}",
                    t as u8, c as u8
                )
            })
            .unwrap_or_else(|| "No bond invariants generator".into());
        Ok(format!(
            "{} --- MorganArguments onlyNonzeroInvariants={} radius={} --- MorganEnvironmentGenerator --- {} --- {}",
            common_arguments_string(
                common.count_simulation,
                common.fp_size,
                common.bits_per_feature,
                common.include_chirality
            ),
            g.only_nonzero_invariants as u8,
            g.radius,
            atoms,
            bonds
        ))
    }
    pub fn to_json(&self) -> Result<String, MorganError> {
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
        // Behavior: full component identities, quoted Boost leaves, source
        // Morgan omission of redundant/bond policies; absent providers omitted.
        // Complexity: O(bounds+query serialization); only explicit serialization
        // invokes the sole Search writer, never fingerprint calculation.
        let state = lock(&self.state)?;
        let mut text = format!(
            "{{\"name\":\"FingerprintGenerator\",\"fingerprintArguments\":{},\"atomEnvironmentGenerator\":{{\"type\":\"MorganEnvGenerator\"}}",
            crate::argument_metadata::morgan_arguments_json(
                state.backend.radius,
                state.backend.only_nonzero_invariants,
                &crate::metadata::common_arguments_json(
                    state.backend.fingerprint_arguments.count_simulation,
                    state.backend.fingerprint_arguments.fp_size,
                    state.backend.fingerprint_arguments.bits_per_feature,
                    state.backend.fingerprint_arguments.include_chirality,
                    &state.backend.fingerprint_arguments.count_bounds
                )
            )
        );
        if let Some(owner) = &state.atom_provider {
            text.push_str(",\"atomInvariantsGenerator\":");
            match owner {
                AtomProvider::Connectivity => text.push_str(&format!(
                    "{{\"type\":\"MorganAtomInvGenerator\",\"includeRingMembership\":\"{}\"}}",
                    state
                        .backend
                        .atom_invariants
                        .as_ref()
                        .expect("connectivity owner")
                        .include_ring_membership()
                )),
                AtomProvider::AtomPair(owner) => text.push_str(&owner.to_json()),
                AtomProvider::Features(owner) => {
                    // RDKit❗✔️: void MorganFeatureAtomInvGenerator::toJSON(
                    // RDKit❗✔️:     boost::property_tree::ptree &pt) const {
                    // RDKit❗✔️:   pt.put("type", "MorganFeatureAtomInvGenerator");
                    // RDKit❗✔️:   if (dp_patterns) {
                    // RDKit❗✔️:     boost::property_tree::ptree patternsNode;
                    // RDKit❗✔️:     for (const auto &pattern : *dp_patterns) {
                    // RDKit❗✔️:       boost::property_tree::ptree patternNode;
                    // RDKit❗✔️:       std::string smarts = MolToSmarts(*pattern);
                    // RDKit❗✔️:       patternNode.put("", smarts);
                    // RDKit❗✔️:       patternsNode.push_back(std::make_pair("", patternNode));
                    // RDKit❗✔️:     }
                    // RDKit❗✔️:     pt.add_child("patternSMARTS", patternsNode);
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   AtomInvariantsGenerator::toJSON(pt);
                    // RDKit❗✔️: }
                    text.push_str("{\"type\":\"MorganFeatureAtomInvGenerator\"");
                    if let Some(patterns) = owner.patterns() {
                        text.push_str(",\"patternSMARTS\":");
                        if patterns.is_empty() {
                            text.push_str("\"\"");
                        } else {
                            text.push('[');
                            for (i, p) in patterns.iter().enumerate() {
                                if i > 0 {
                                    text.push(',');
                                }
                                let smarts =
                                    query_graph_to_smarts(p, &SmartsWriteParams::default())?;
                                text.push_str(&Value::String(smarts).to_string());
                            }
                            text.push(']');
                        }
                    }
                    text.push('}');
                }
            }
        }
        if let Some(owner) = &state.backend.bond_invariants {
            let (t, c) = owner.configuration();
            text.push_str(&format!(",\"bondInvariantsGenerator\":{{\"type\":\"MorganBondInvGenerator\",\"useBondTypes\":\"{}\",\"useChirality\":\"{}\"}}",t,c));
        }
        text.push('}');
        Ok(text)
    }
    pub fn from_json(json: &str) -> Result<Self, MorganError> {
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
        // Behavior: default source args updated without repeating constructor
        // checks; missing invariant providers remain absent, not defaults.
        // Known unmodeled RDKitFP providers report their independent capability.
        // Complexity: one JSON parse and one optional query preparation pass.
        let value = parse_object(json)?;
        let child = |name: &str| {
            value.get(name).filter(|v| v.is_object()).ok_or_else(|| {
                FingerprintJsonError::Invalid(format!("missing {name} node in JSON"))
            })
        };
        fn typ<'a>(node: &'a SourceNode, category: &str) -> Result<&'a str, FingerprintJsonError> {
            node.get("type")
                .and_then(SourceNode::as_str)
                .ok_or_else(|| {
                    FingerprintJsonError::Invalid(format!("{category} type not specified in JSON"))
                })
        }
        let args = child("fingerprintArguments")?;
        if typ(args, "FingerprintArguments")? != "MorganArguments" {
            return Err(FingerprintJsonError::Invalid(
                "JSON does not describe a Morgan fingerprint generator".into(),
            )
            .into());
        }
        let env = child("atomEnvironmentGenerator")?;
        let env_type = typ(env, "AtomEnvironmentGenerator")?;
        if env_type != "MorganEnvGenerator" {
            return Err(if matches!(
                env_type,
                "RDKitFPEnvGenerator" | "AtomPairEnvGenerator" | "TopologicalTorsionEnvGenerator"
            ) {
                FingerprintJsonError::UnsupportedComponent {
                    component: "atomEnvironmentGenerator",
                    source_type: env_type.into(),
                }
            } else {
                FingerprintJsonError::Invalid(format!("Unknown AtomEnvGenerator type: {env_type}"))
            }
            .into());
        }
        let params = MorganParams::default().with_json_value(args)?;
        let mut backend = MorganGenerator {
            radius: params.radius,
            only_nonzero_invariants: params.only_nonzero_invariants,
            include_redundant_environments: params.include_redundant_environments,
            fingerprint_arguments: FingerprintArguments {
                count_simulation: params.count_simulation,
                include_chirality: params.include_chirality,
                count_bounds: params.count_bounds,
                fp_size: params.fp_size,
                bits_per_feature: params.bits_per_feature,
            },
            atom_invariants: None,
            atom_provider_present: false,
            bond_invariants: None,
        };
        let atom_provider = if let Some(node) = value.get("atomInvariantsGenerator") {
            let atom_type = typ(node, "AtomInvariantsGenerator")?;
            backend.atom_provider_present = true;
            Some(match atom_type {
                "MorganAtomInvGenerator" => {
                    let ring = bool_or(node, "includeRingMembership", true);
                    backend.atom_invariants = Some(MorganAtomInvGenerator::new(ring));
                    AtomProvider::Connectivity
                }
                "AtomPairAtomInvGenerator" => {
                    let mut owner = AtomPairAtomInvariantsGenerator::default();
                    owner.from_json_value(node)?;
                    AtomProvider::AtomPair(owner)
                }
                "MorganFeatureAtomInvGenerator" => {
                    // RDKit❗✔️: void MorganFeatureAtomInvGenerator::fromJSON(
                    // RDKit❗✔️:     const boost::property_tree::ptree &pt) {
                    // RDKit❗✔️:   if (pt.get_child_optional("patternSMARTS")) {
                    // RDKit❗✔️:     const auto &patternsNode = pt.get_child("patternSMARTS");
                    // RDKit❗✔️:     cleanUpPatterns();
                    // RDKit❗✔️:     dp_patterns = new std::vector<const ROMol *>();
                    // RDKit❗✔️:     for (const auto &patternNode : patternsNode) {
                    // RDKit❗✔️:       std::string smarts = patternNode.second.get_value<std::string>();
                    // RDKit❗✔️:       ROMol *patternMol = SmartsToMol(smarts);
                    // RDKit❗✔️:       if (patternMol) {
                    // RDKit❗✔️:         dp_patterns->push_back(patternMol);
                    // RDKit❗✔️:       }
                    // RDKit❗✔️:     }
                    // RDKit❗✔️:   }
                    // RDKit❗✔️:   AtomInvariantsGenerator::fromJSON(pt);
                    // RDKit❗✔️: }
                    let patterns = if let Some(field) = node.get("patternSMARTS") {
                        let mut queries = Vec::new();
                        {
                            for pattern in field.children() {
                                let smarts = pattern.as_str().ok_or_else(|| {
                                    FingerprintJsonError::Invalid(
                                        "patternSMARTS entries must be strings".into(),
                                    )
                                })?;
                                // Source nullptr boundary: exact typed syntax is resolved by
                                // the sole SEARCH parser. Generic Parse still mixes unresolved
                                // syntax with internal graph errors and remains propagated.
                                // RDKit❗✔️: std::unique_ptr<RWMol> toMol(const std::string &inp,
                                // RDKit❗✔️:                              int func(const std::string &,
                                // RDKit❗✔️:                                       std::vector<RDKit::RWMol *> &),
                                // RDKit❗✔️:                              const std::string &origInp) {
                                // RDKit❗✔️:   // empty strings produce empty molecules:
                                // RDKit❗✔️:   if (inp.empty()) {
                                // RDKit❗✔️:     return std::make_unique<RWMol>();
                                // RDKit❗✔️:   }
                                // RDKit❗✔️:   std::unique_ptr<RWMol> res;
                                // RDKit❗✔️:   std::vector<RDKit::RWMol *> molVect;
                                // RDKit❗✔️:   try {
                                // RDKit❗✔️:     func(inp, molVect);
                                // RDKit❗✔️:     if (!molVect.empty()) {
                                // RDKit❗✔️:       res.reset(molVect[0]);
                                // RDKit❗✔️:       SmilesParseOps::CloseMolRings(res.get(), false);
                                // RDKit❗✔️:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
                                // RDKit❗✔️:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
                                // RDKit❗✔️:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
                                // RDKit❗✔️:       // No sense leaving this bookmark intact:
                                // RDKit❗✔️:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
                                // RDKit❗✔️:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
                                // RDKit❗✔️:       }
                                // RDKit❗✔️:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
                                // RDKit❗✔️:                              // occur last in this if.
                                // RDKit❗✔️:     }
                                // RDKit❗✔️:   } catch (SmilesParseException &e) {
                                // RDKit❗✔️:     std::string nm = "SMILES";
                                // RDKit❗✔️:     if (func == smarts_parse) {
                                // RDKit❗✔️:       nm = "SMARTS";
                                // RDKit❗✔️:     }
                                // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
                                // RDKit❗✔️:                           << " for input: '" << origInp << "'" << std::endl;
                                // RDKit❗✔️:
                                // RDKit❗✔️:     // reset res so that we return a nullptr. We don't want to reset(),
                                // RDKit❗✔️:     // because that would delete the mol and leak any unmatched
                                // RDKit❗✔️:     // ring closure bonds. These will be cleaned up in the loop below.
                                // RDKit❗✔️:     res.release();
                                // RDKit❗✔️:   }
                                // RDKit❗✔️:   for (auto *molPtr : molVect) {
                                // RDKit❗✔️:     if (molPtr) {
                                // RDKit❗✔️:       // Clean-up the bond bookmarks when not calling CloseMolRings
                                // RDKit❗✔️:       SmilesParseOps::CleanupAfterParseError(molPtr);
                                // RDKit❗✔️:       delete molPtr;
                                // RDKit❗✔️:     }
                                // RDKit❗✔️:   }
                                // RDKit❗✔️:
                                // RDKit❗✔️:   return res;
                                // RDKit❗✔️: }
                                match parse_smarts(smarts, &SmartsParseParams::default()) {
                                    Ok(query) => queries.push(query),
                                    Err(
                                        SmartsParseError::UnclosedBracket(_)
                                        | SmartsParseError::UnclosedParenthesis(_)
                                        | SmartsParseError::UnexpectedEnd(_)
                                        | SmartsParseError::UnexpectedCharacter { .. }
                                        | SmartsParseError::InvalidAtomPrimitive { .. }
                                        | SmartsParseError::UnbalancedRingClosure(_)
                                        | SmartsParseError::SelfRingClosure { .. }
                                        | SmartsParseError::DuplicateRingBond { .. },
                                    ) => {}
                                    Err(error) => return Err(error.into()),
                                }
                            }
                        }
                        Some(queries)
                    } else {
                        None
                    };
                    AtomProvider::Features(MorganFeatureAtomInvGenerator::from_owned_patterns(
                        patterns,
                    ))
                }
                "RDKitFPAtomInvGenerator" => {
                    return Err(FingerprintJsonError::UnsupportedComponent {
                        component: "atomInvariantsGenerator",
                        source_type: atom_type.into(),
                    }
                    .into());
                }
                _ => {
                    return Err(FingerprintJsonError::Invalid(format!(
                        "Unknown AtomInvariantsGenerator type: {atom_type}"
                    ))
                    .into());
                }
            })
        } else {
            None
        };
        if let Some(node) = value.get("bondInvariantsGenerator") {
            let bond_type = typ(node, "BondInvariantsGenerator")?;
            if bond_type != "MorganBondInvGenerator" {
                return Err(FingerprintJsonError::Invalid(format!(
                    "Unknown BondInvariantsGenerator type: {bond_type}"
                ))
                .into());
            }
            let t = bool_or(node, "useBondTypes", true);
            let c = bool_or(node, "useChirality", false);
            backend.bond_invariants = Some(MorganBondInvGenerator::new(t, c));
        }
        Ok(Self {
            state: Arc::new(Mutex::new(State {
                backend,
                atom_provider,
                initial_ring_membership: true,
                initial_use_bond_types: true,
            })),
        })
    }
}
macro_rules! calculation {
    ($name:ident,$state_name:ident,$source:ident,$result:ty,$bulk:ident) => {
        impl MorganOperator {
            pub fn $name(
                &self,
                input: &MorganPreparedInput<'_>,
                call: &MorganCall<'_>,
                output: Option<&mut FingerprintAdditionalOutput>,
            ) -> Result<$result, MorganError> {
                let state = lock(&self.state)?;
                Self::$state_name(&state, input, call, output)
            }
            fn $state_name(
                state: &State,
                input: &MorganPreparedInput<'_>,
                call: &MorganCall<'_>,
                output: Option<&mut FingerprintAdditionalOutput>,
            ) -> Result<$result, MorganError> {
                // One source dispatch to the persistent configured backend;
                // common preparation, custom precedence and projection remain
                // in the existing source owner, without per-call factory work.
                let args = call.fingerprint_func_arguments();
                $source(
                    input.topology,
                    input.properties,
                    input.valence,
                    input.rings,
                    &state.backend,
                    &args,
                    output,
                    |_t, _p, _v, _r, _g| state.atoms(input),
                )
            }
            pub fn $bulk(
                &self,
                inputs: &[Option<MorganPreparedInput<'_>>],
                num_threads: i32,
            ) -> Result<Vec<Option<$result>>, MorganError> {
                let state = lock(&self.state)?;
                let workers = crate::fingerprint_bulk::workers(num_threads)?;
                let call = MorganCall::default();
                crate::fingerprint_bulk::ordered_bulk(inputs, workers, |input| {
                    Self::$state_name(&state, input, &call, None)
                })
            }
        }
    };
}
calculation!(
    sparse_count,
    sparse_count_state,
    get_sparse_count_fingerprint_with_atom_invariants,
    crate::SparseCountFingerprint,
    sparse_counts
);
calculation!(
    sparse_bits,
    sparse_bits_state,
    get_sparse_fingerprint_with_atom_invariants,
    crate::SparseBitFingerprint,
    sparse_fingerprints
);
calculation!(
    count,
    count_state,
    get_count_fingerprint_with_atom_invariants,
    SparseCountFingerprint32,
    counts
);
calculation!(
    bits,
    bits_state,
    get_fingerprint_with_atom_invariants,
    Fingerprint,
    fingerprints
);
macro_rules! setting {($get:ident,$set:ident,$ty:ty,$($field:ident).+) => {
    pub fn $get(&self)->Result<$ty,MorganError>{Ok(lock(&self.state)?.backend.$($field).+)}
    pub fn $set(&self,value:$ty)->Result<(),MorganError>{lock(&self.state)?.backend.$($field).+ = value;Ok(())}
}}
impl MorganSettings {
    pub fn snapshot(&self) -> Result<MorganParams, MorganError> {
        Ok(lock(&self.state)?.params())
    }
    setting!(radius, set_radius, u32, radius);
    setting!(
        only_nonzero_invariants,
        set_only_nonzero_invariants,
        bool,
        only_nonzero_invariants
    );
    setting!(
        include_redundant_environments,
        set_include_redundant_environments,
        bool,
        include_redundant_environments
    );
    setting!(
        count_simulation,
        set_count_simulation,
        bool,
        fingerprint_arguments.count_simulation
    );
    setting!(
        include_chirality,
        set_include_chirality,
        bool,
        fingerprint_arguments.include_chirality
    );
    setting!(fp_size, set_fp_size, u32, fingerprint_arguments.fp_size);
    setting!(
        bits_per_feature,
        set_bits_per_feature,
        u32,
        fingerprint_arguments.bits_per_feature
    );
    pub fn count_bounds(&self) -> Result<Vec<u32>, MorganError> {
        Ok(lock(&self.state)?
            .backend
            .fingerprint_arguments
            .count_bounds
            .clone())
    }
    pub fn set_count_bounds(&self, bounds: Vec<u32>) -> Result<(), MorganError> {
        lock(&self.state)?
            .backend
            .fingerprint_arguments
            .count_bounds = bounds;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestMolecule;
    use std::collections::BTreeMap;
    fn input(m: &TestMolecule) -> MorganPreparedInput<'_> {
        MorganPreparedInput {
            topology: &m.topology,
            properties: &m.properties,
            coordinates: &m.coordinates,
            valence: &m.valence,
            rings: &m.rings,
        }
    }
    fn generator() -> MorganOperator {
        MorganOperator::new(&MorganParams::default(), None, None).unwrap()
    }
    #[test]
    fn source_self_ring_and_duplicate_bond_null_patterns_preserve_valid_query_order() {
        let mut value: Value = serde_json::from_str(&generator().to_json().unwrap()).unwrap();
        value["fingerprintArguments"]["radius"] = Value::String("0".into());
        let molecule = TestMolecule::from_smiles("CCO").unwrap();
        for invalid in ["[C]11", "[C]1[C]1"] {
            value["atomInvariantsGenerator"] = serde_json::json!({
                "type": "MorganFeatureAtomInvGenerator", "patternSMARTS": ["[C]", invalid, "[O]"]
            });
            let restored = MorganOperator::from_json(&value.to_string()).unwrap();
            let json: Value = serde_json::from_str(&restored.to_json().unwrap()).unwrap();
            assert_eq!(
                json["atomInvariantsGenerator"]["patternSMARTS"],
                serde_json::json!(["C", "O"])
            );
            assert_eq!(
                restored
                    .sparse_count(&input(&molecule), &MorganCall::default(), None)
                    .unwrap()
                    .nonzero_elements(),
                &BTreeMap::from([(1, 2), (2, 1)])
            );
        }
        value["atomInvariantsGenerator"]["patternSMARTS"] =
            serde_json::json!(["[C]11", "[C]1[C]1"]);
        let restored = MorganOperator::from_json(&value.to_string()).unwrap();
        assert_eq!(
            restored
                .sparse_count(&input(&molecule), &MorganCall::default(), None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(0, 3)])
        );
    }

    #[test]
    fn generic_source_parser_category_is_propagated_without_string_classification() {
        let mut value: Value = serde_json::from_str(&generator().to_json().unwrap()).unwrap();
        value["atomInvariantsGenerator"] = serde_json::json!({
            "type": "MorganFeatureAtomInvGenerator", "patternSMARTS": ["[C]", "C1CC", "[O]"]
        });
        let error = MorganOperator::from_json(&value.to_string()).unwrap_err();
        assert!(matches!(
            &error,
            MorganError::SmartsParse(SmartsParseError::Parse(_))
        ));
        let source = std::error::Error::source(&error).unwrap();
        assert!(matches!(
            source.downcast_ref::<SmartsParseError>(),
            Some(SmartsParseError::Parse(_))
        ));
        // This preserves ROOT's frozen generic-error boundary, not an assertion
        // that the unresolved unclosed-ring/null syntax distinction is complete.
    }

    #[test]
    fn fixed_complete_source_metadata_and_cco_raw_count_output() {
        let g = generator();
        assert_eq!(
            g.info_string().unwrap(),
            "Common arguments : countSimulation=0 fpSize=2048 bitsPerFeature=1 includeChirality=0 --- MorganArguments onlyNonzeroInvariants=0 radius=3 --- MorganEnvironmentGenerator --- MorganInvariantGenerator includeRingMembership=1 --- MorganInvariantGenerator useBondTypes=1 useChirality=0"
        );
        let value: Value = serde_json::from_str(&g.to_json().unwrap()).unwrap();
        assert_eq!(
            value,
            serde_json::json!({"name":"FingerprintGenerator","fingerprintArguments":{"type":"MorganArguments","onlyNonzeroInvariants":"false","radius":"3","countSimulation":"false","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]},"atomEnvironmentGenerator":{"type":"MorganEnvGenerator"},"atomInvariantsGenerator":{"type":"MorganAtomInvGenerator","includeRingMembership":"true"},"bondInvariantsGenerator":{"type":"MorganBondInvGenerator","useBondTypes":"true","useChirality":"false"}})
        );
        let m = TestMolecule::from_smiles("CCO").unwrap();
        assert_eq!(
            g.sparse_count(&input(&m), &MorganCall::default(), None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([
                (864662311, 1),
                (1535166686, 1),
                (2245384272, 1),
                (2246728737, 1),
                (3542456614, 1),
                (4018048386, 1)
            ])
        );
    }
    #[test]
    fn shared_live_settings_survive_operator_alias_lifetime_and_preserve_captured_provider() {
        let g = generator();
        let alias = g.clone();
        let settings = g.settings();
        let other = alias.settings();
        settings.set_radius(0).unwrap();
        settings.set_include_chirality(true).unwrap();
        settings.set_count_bounds(vec![1, 3, 7]).unwrap();
        assert_eq!(other.radius().unwrap(), 0);
        assert!(other.include_chirality().unwrap());
        assert_eq!(other.count_bounds().unwrap(), [1, 3, 7]);
        assert!(alias.info_string().unwrap().contains("includeChirality=1"));
        assert!(
            alias
                .info_string()
                .unwrap()
                .contains("useBondTypes=1 useChirality=0")
        );
        let m = TestMolecule::from_smiles("CCO").unwrap();
        assert_eq!(
            alias
                .sparse_count(&input(&m), &MorganCall::default(), None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(864662311, 1), (2245384272, 1), (2246728737, 1)])
        );
        drop(g);
        drop(alias);
        settings.set_radius(2).unwrap();
        assert_eq!(other.radius().unwrap(), 2);
        let mut snapshot = other.snapshot().unwrap();
        snapshot.radius = 77;
        snapshot.count_bounds.push(99);
        assert_eq!(settings.radius().unwrap(), 2);
        assert_eq!(settings.count_bounds().unwrap(), [1, 3, 7]);
    }
    #[test]
    fn explicit_provider_precedence_ring_configuration_and_live_chirality_capture() {
        let params = MorganParams {
            include_ring_membership: false,
            use_bond_types: true,
            include_chirality: false,
            ..Default::default()
        };
        let default = MorganOperator::new(&params, None, None).unwrap();
        assert!(
            default
                .info_string()
                .unwrap()
                .contains("includeRingMembership=0")
        );
        let explicit = MorganOperator::new(
            &params,
            Some(MorganAtomProvider::Connectivity {
                include_ring_membership: true,
            }),
            Some(MorganBondProvider {
                use_bond_types: false,
                include_chirality: true,
            }),
        )
        .unwrap();
        assert!(
            explicit
                .info_string()
                .unwrap()
                .contains("includeRingMembership=1")
        );
        assert!(
            explicit
                .info_string()
                .unwrap()
                .contains("useBondTypes=0 useChirality=1")
        );
        explicit.settings().set_include_chirality(true).unwrap();
        explicit.settings().set_include_chirality(false).unwrap();
        assert!(
            explicit
                .info_string()
                .unwrap()
                .contains("useBondTypes=0 useChirality=1")
        );
        let ring = TestMolecule::from_smiles("C1CCCCC1").unwrap();
        assert_ne!(
            default
                .sparse_count(&input(&ring), &MorganCall::default(), None)
                .unwrap(),
            explicit
                .sparse_count(&input(&ring), &MorganCall::default(), None)
                .unwrap()
        );
    }
    #[test]
    fn four_outputs_and_full_additional_output_reuse_existing_source_pipeline() {
        for text in ["CCO", "C1CCCCC1", "C[C@H](F)Cl"] {
            let m = TestMolecule::from_smiles(text).unwrap();
            let p = MorganParams::default();
            let g = MorganOperator::new(&p, None, None).unwrap();
            let call = MorganCall::default();
            let i = input(&m);
            macro_rules! compare {
                ($method:ident,$original:ident) => {{
                    let mut actual = FingerprintAdditionalOutput::new();
                    let mut expected = FingerprintAdditionalOutput::new();
                    for out in [&mut actual, &mut expected] {
                        out.allocate_atom_counts();
                        out.allocate_atom_to_bits();
                        out.allocate_bit_info_map();
                        out.allocate_bit_paths();
                        out.allocate_atoms_per_bit();
                    }
                    assert_eq!(
                        g.$method(&i, &call, Some(&mut actual)).unwrap(),
                        $original(
                            &i,
                            &p,
                            &call,
                            MorganAtomInvariants::Connectivity,
                            Some(&mut expected)
                        )
                        .unwrap()
                    );
                    assert_eq!(actual, expected);
                }};
            }
            compare!(sparse_count, morgan_sparse_count);
            compare!(sparse_bits, morgan_sparse_bits);
            compare!(count, morgan_count);
            compare!(bits, morgan_bits);
        }
    }
    #[test]
    fn source_feature_patterns_present_empty_and_custom_invariants_have_distinct_semantics() {
        let p = MorganParams {
            radius: 0,
            ..Default::default()
        };
        let m = TestMolecule::from_smiles("CCO").unwrap();
        let call = MorganCall::default();
        let empty = MorganOperator::new(
            &p,
            Some(MorganAtomProvider::Features {
                patterns: Some(vec![]),
            }),
            None,
        )
        .unwrap();
        assert_eq!(
            empty
                .sparse_count(&input(&m), &call, None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(0, 3)])
        );
        let v: Value = serde_json::from_str(&empty.to_json().unwrap()).unwrap();
        assert_eq!(v["atomInvariantsGenerator"]["patternSMARTS"], "");
        let c = parse_smarts("[C]", &SmartsParseParams::default()).unwrap();
        let o = parse_smarts("[O]", &SmartsParseParams::default()).unwrap();
        let configured = MorganOperator::new(
            &p,
            Some(MorganAtomProvider::Features {
                patterns: Some(vec![c, o]),
            }),
            None,
        )
        .unwrap();
        assert_eq!(
            configured
                .sparse_count(&input(&m), &call, None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(1, 2), (2, 1)])
        );
        let restored = MorganOperator::from_json(&configured.to_json().unwrap()).unwrap();
        assert_eq!(
            restored.sparse_count(&input(&m), &call, None).unwrap(),
            configured.sparse_count(&input(&m), &call, None).unwrap()
        );
        let override_call = MorganCall {
            custom_atom_invariants: Some(&[11, 12, 13]),
            ..Default::default()
        };
        assert_eq!(
            configured
                .sparse_count(&input(&m), &override_call, None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(11, 1), (12, 1), (13, 1)])
        );
        assert!(matches!(
            configured.sparse_count(
                &input(&m),
                &MorganCall {
                    custom_atom_invariants: Some(&[]),
                    ..Default::default()
                },
                None
            ),
            Err(MorganError::Fingerprint(
                FingerprintError::PreconditionViolation {
                    what: "bad atom invariants size"
                }
            ))
        ));
    }
    #[test]
    fn source_feature_json_invalid_pattern_null_branch_and_provider_omission() {
        let g = MorganOperator::new(
            &MorganParams {
                radius: 0,
                ..Default::default()
            },
            Some(MorganAtomProvider::Features { patterns: None }),
            None,
        )
        .unwrap();
        let mut v: Value = serde_json::from_str(&g.to_json().unwrap()).unwrap();
        assert!(v["atomInvariantsGenerator"].get("patternSMARTS").is_none());
        v["atomInvariantsGenerator"]["patternSMARTS"] = serde_json::json!(["[C]", "[", "[O]"]);
        let restored = MorganOperator::from_json(&v.to_string()).unwrap();
        let m = TestMolecule::from_smiles("CCO").unwrap();
        assert_eq!(
            restored
                .sparse_count(&input(&m), &MorganCall::default(), None)
                .unwrap()
                .nonzero_elements(),
            &BTreeMap::from([(1, 2), (2, 1)])
        );
        assert_eq!(serde_json::from_str::<Value>(&restored.to_json().unwrap()).unwrap()["atomInvariantsGenerator"]["patternSMARTS"].as_array().unwrap().len(),2);
    }
    #[test]
    fn absent_json_providers_preserve_null_preconditions_even_for_empty_molecule() {
        let g = generator();
        let mut value: Value = serde_json::from_str(&g.to_json().unwrap()).unwrap();
        value
            .as_object_mut()
            .unwrap()
            .remove("atomInvariantsGenerator");
        value
            .as_object_mut()
            .unwrap()
            .remove("bondInvariantsGenerator");
        let restored = MorganOperator::from_json(&value.to_string()).unwrap();
        let empty = TestMolecule::new();
        let i = input(&empty);
        assert!(
            restored
                .info_string()
                .unwrap()
                .contains("No atom invariants generator --- No bond invariants generator")
        );
        assert!(matches!(
            restored.sparse_count(&i, &MorganCall::default(), None),
            Err(MorganError::Fingerprint(
                FingerprintError::PreconditionViolation {
                    what: "bad atom invariants size"
                }
            ))
        ));
        assert!(matches!(
            restored.sparse_count(
                &i,
                &MorganCall {
                    custom_atom_invariants: Some(&[]),
                    ..Default::default()
                },
                None
            ),
            Err(MorganError::Fingerprint(
                FingerprintError::PreconditionViolation {
                    what: "bad bond invariants size"
                }
            ))
        ));
        assert!(
            restored
                .sparse_count(
                    &i,
                    &MorganCall {
                        custom_atom_invariants: Some(&[]),
                        custom_bond_invariants: Some(&[]),
                        ..Default::default()
                    },
                    None
                )
                .unwrap()
                .nonzero_elements()
                .is_empty()
        );
        let m = TestMolecule::from_smiles("CCO").unwrap();
        for threads in [1, 3] {
            assert!(matches!(
                restored.sparse_counts(&[None, Some(input(&m))], threads),
                Err(MorganError::Fingerprint(
                    FingerprintError::PreconditionViolation {
                        what: "bad atom invariants size"
                    }
                ))
            ));
        }
    }
    #[test]
    fn source_json_roundtrip_and_omitted_arguments_do_not_refresh_provider_or_validate_live_settings()
     {
        let g = generator();
        g.settings().set_include_chirality(true).unwrap();
        g.settings()
            .set_include_redundant_environments(true)
            .unwrap();
        g.settings().set_count_bounds(vec![]).unwrap();
        g.settings().set_bits_per_feature(0).unwrap();
        let restored = MorganOperator::from_json(&g.to_json().unwrap()).unwrap();
        assert!(restored.settings().include_chirality().unwrap());
        assert!(
            !restored
                .settings()
                .include_redundant_environments()
                .unwrap()
        );
        assert_eq!(restored.settings().bits_per_feature().unwrap(), 0);
        assert!(restored.settings().count_bounds().unwrap().is_empty());
        assert!(restored.info_string().unwrap().contains("useChirality=0"));
        let m = TestMolecule::from_smiles("CCO").unwrap();
        assert_eq!(
            restored
                .sparse_count(&input(&m), &MorganCall::default(), None)
                .unwrap()
                .nonzero_elements()
                .len(),
            6
        );
    }
    #[test]
    fn source_bulk_four_outputs_order_missing_and_empty_slots() {
        let g = generator();
        let a = TestMolecule::from_smiles("CCO").unwrap();
        let b = TestMolecule::from_smiles("C").unwrap();
        let e = TestMolecule::new();
        let rows = [
            Some(input(&a)),
            None,
            Some(input(&b)),
            Some(input(&e)),
            None,
        ];
        for threads in [1, 2, 3, 7] {
            macro_rules! compare {
                ($bulk:ident,$single:ident) => {{
                    let expected = rows
                        .iter()
                        .map(|r| {
                            r.as_ref()
                                .map(|i| g.$single(i, &MorganCall::default(), None).unwrap())
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(g.$bulk(&rows, threads).unwrap(), expected);
                    assert!(g.$bulk(&[], threads).unwrap().is_empty());
                }};
            }
            compare!(sparse_counts, sparse_count);
            compare!(sparse_fingerprints, sparse_bits);
            compare!(counts, count);
            compare!(fingerprints, bits);
        }
    }
    #[test]
    fn known_independent_unmodeled_component_and_unknown_type_errors_remain_distinct() {
        let g = generator();
        let mut value: Value = serde_json::from_str(&g.to_json().unwrap()).unwrap();
        value["atomInvariantsGenerator"]["type"] = Value::String("RDKitFPAtomInvGenerator".into());
        assert!(matches!(
            MorganOperator::from_json(&value.to_string()),
            Err(MorganError::Json(
                FingerprintJsonError::UnsupportedComponent {
                    component: "atomInvariantsGenerator",
                    ..
                }
            ))
        ));
        value["atomInvariantsGenerator"]["type"] = Value::String("invented".into());
        assert!(matches!(
            MorganOperator::from_json(&value.to_string()),
            Err(MorganError::Json(FingerprintJsonError::Invalid(_)))
        ));
        assert!(matches!(
            MorganOperator::from_json("not json"),
            Err(MorganError::Json(FingerprintJsonError::Parse(_)))
        ));
    }
}
