//! Canonical cross-language naming registry.

use cosmolkit_macros::binding_contract;

binding_contract! {
    pub static BINDING_CONTRACT = [
        {
            semantic_id: "types.Element", item: type, owner: type_,
            rust: crate::Element, python: "Element", javascript: "Element",
            feature: "metadata", status: experimental, role: value,
        },
        {
            semantic_id: "Element.from_atomic_number", item: callable, owner: type_,
            rust: crate::Element::from_atomic_number, python: "from_atomic_number", javascript: "fromAtomicNumber",
            feature: "metadata", status: experimental, kind: static_,
            parameters: [{ name: atomic_number, type: u8, default: required }],
            output: Option<crate::Element>, error: none,
            state: value_returning, operation: none,
            signature: fn(u8) -> Option<crate::Element>,
        },
        {
            semantic_id: "Element.from_symbol", item: callable, owner: type_,
            rust: crate::Element::from_symbol, python: "from_symbol", javascript: "fromSymbol",
            feature: "metadata", status: experimental, kind: static_,
            parameters: [{ name: symbol, type: &str, default: required }],
            output: Option<crate::Element>, error: none,
            state: value_returning, operation: none,
            signature: for<'a> fn(&'a str) -> Option<crate::Element>,
        },
        {
            semantic_id: "Element.atomic_number", item: callable, owner: type_,
            rust: crate::Element::atomic_number, python: "atomic_number", javascript: "atomicNumber",
            feature: "metadata", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: u8, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::Element) -> u8,
        },
        {
            semantic_id: "Element.symbol", item: callable, owner: type_,
            rust: crate::Element::symbol, python: "symbol", javascript: "symbol",
            feature: "metadata", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: &'static str, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::Element) -> &'static str,
        },
        {
            semantic_id: "types.ElementInfo", item: type, owner: type_,
            rust: crate::ElementInfo, python: "ElementInfo", javascript: "ElementInfo",
            feature: "metadata", status: experimental, role: result,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "module.element_info", item: callable, owner: module,
            rust: crate::element_info, python: "element_info", javascript: "elementInfo",
            feature: "cap-valence", status: experimental, kind: module,
            parameters: [{ name: element, type: crate::Element, default: required }],
            output: crate::ElementInfo, error: none,
            state: read_only, operation: none,
            signature: fn(crate::Element) -> crate::ElementInfo,
        },
        {
            semantic_id: "types.QueryGraph", item: type, owner: type_,
            rust: crate::QueryGraph, python: "QueryGraph", javascript: "QueryGraph",
            feature: "metadata", status: experimental, role: value,
        },
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"errors.TopologicalTorsionPathScoreError",item:type,owner:type_,rust:crate::TopologicalTorsionPathScoreError,python:"TopologicalTorsionPathScoreError",javascript:"TopologicalTorsionPathScoreError",feature:"cap-fingerprints",status:experimental,role:error,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.topological_torsion_path_score",item:callable,owner:molecule,rust:crate::Molecule::topological_torsion_path_score,python:"topological_torsion_path_score",javascript:"topologicalTorsionPathScore",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:path,type:&[usize],default:required},{name:size,type:usize,default:required},{name:atom_codes,type:Option<&[u32]>,default:none}],output:u64,error:crate::TopologicalTorsionPathScoreError,state:read_only,operation:none,signature:for<'a,'b,'c> fn(&'a crate::Molecule,&'b [usize],usize,Option<&'c [u32]>)->Result<u64,crate::TopologicalTorsionPathScoreError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"explain_path_score",item:callable,owner:module,rust:crate::explain_path_score,python:"explain_path_score",javascript:"explainPathScore",feature:"cap-fingerprints",status:experimental,kind:module,parameters:[{name:score,type:u64,default:required},{name:size,type:usize,default:integer(4)}],output:Vec<(&'static str,u32,u32)>,error:none,state:read_only,operation:none,signature:fn(u64,usize)->Vec<(&'static str,u32,u32)>,},
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomPairsParameters", item: type, owner: type_, rust: crate::AtomPairsParameters, python: "AtomPairsParameters", javascript: "AtomPairsParameters", feature: "cap-fingerprints", status: experimental, role: value, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.version", item: callable, owner: type_, rust: crate::AtomPairsParameters::version, python: "version", javascript: "version", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: &'static str, error: none, state: read_only, operation: none, signature: fn()->&'static str, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_type_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_type_bits, python: "num_type_bits", javascript: "numTypeBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_pi_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_pi_bits, python: "num_pi_bits", javascript: "numPiBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_branch_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_branch_bits, python: "num_branch_bits", javascript: "numBranchBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_chiral_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_chiral_bits, python: "num_chiral_bits", javascript: "numChiralBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.code_size", item: callable, owner: type_, rust: crate::AtomPairsParameters::code_size, python: "code_size", javascript: "codeSize", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_path_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_path_bits, python: "num_path_bits", javascript: "numPathBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.max_path_length", item: callable, owner: type_, rust: crate::AtomPairsParameters::max_path_length, python: "max_path_length", javascript: "maxPathLength", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.num_atom_pair_fingerprint_bits", item: callable, owner: type_, rust: crate::AtomPairsParameters::num_atom_pair_fingerprint_bits, python: "num_atom_pair_fingerprint_bits", javascript: "numAtomPairFingerprintBits", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn()->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairsParameters.atom_types", item: callable, owner: type_, rust: crate::AtomPairsParameters::atom_types, python: "atom_types", javascript: "atomTypes", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [], output: Vec<u32>, error: none, state: read_only, operation: none, signature: fn()->Vec<u32>, },






































        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.MorganAtomInvariantsGenerator",item:type,owner:type_,rust:crate::MorganAtomInvariantsGenerator,python:"MorganAtomInvariantsGenerator",javascript:"MorganAtomInvariantsGenerator",feature:"cap-fingerprints",status:experimental,role:parameter,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.MorganBondInvariantsGenerator",item:type,owner:type_,rust:crate::MorganBondInvariantsGenerator,python:"MorganBondInvariantsGenerator",javascript:"MorganBondInvariantsGenerator",feature:"cap-fingerprints",status:experimental,role:parameter,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.MorganFingerprintGenerator",item:type,owner:type_,rust:crate::MorganFingerprintGenerator,python:"MorganFingerprintGenerator",javascript:"MorganFingerprintGenerator",feature:"cap-fingerprints",status:experimental,role:value,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.MorganSettings",item:type,owner:type_,rust:crate::MorganSettings,python:"MorganSettings",javascript:"MorganSettings",feature:"cap-fingerprints",status:experimental,role:value,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.MorganCallParams",item:type,owner:type_,rust:crate::MorganCallParams,python:"MorganCallParams",javascript:"MorganCallParams",feature:"cap-fingerprints",status:experimental,role:parameter,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganAtomInvariantsGenerator.connectivity",item:callable,owner:type_,rust:crate::MorganAtomInvariantsGenerator::connectivity,python:"connectivity",javascript:"connectivity",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:include_ring_membership,type:bool,default:"true"}],output:crate::MorganAtomInvariantsGenerator,error:none,state:value_returning,operation:none,signature:fn(bool)->crate::MorganAtomInvariantsGenerator,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganAtomInvariantsGenerator.features",item:callable,owner:type_,rust:crate::MorganAtomInvariantsGenerator::features,python:"features",javascript:"features",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:patterns,type:Option<Vec<crate::QueryGraph>>,default:none}],output:crate::MorganAtomInvariantsGenerator,error:none,state:value_returning,operation:none,signature:fn(Option<Vec<crate::QueryGraph>>)->crate::MorganAtomInvariantsGenerator,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganAtomInvariantsGenerator.atom_pair",item:callable,owner:type_,rust:crate::MorganAtomInvariantsGenerator::atom_pair,python:"atom_pair",javascript:"atomPair",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:generator,type:crate::AtomPairAtomInvariantsGenerator,default:required}],output:crate::MorganAtomInvariantsGenerator,error:none,state:value_returning,operation:none,signature:fn(crate::AtomPairAtomInvariantsGenerator)->crate::MorganAtomInvariantsGenerator,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganBondInvariantsGenerator.new",item:callable,owner:type_,rust:crate::MorganBondInvariantsGenerator::new,python:"new",javascript:"new",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:use_bond_types,type:bool,default:"true"},{name:include_chirality,type:bool,default:"false"}],output:crate::MorganBondInvariantsGenerator,error:none,state:value_returning,operation:none,signature:fn(bool,bool)->crate::MorganBondInvariantsGenerator,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganBondInvariantsGenerator.use_bond_types",item:callable,owner:type_,rust:crate::MorganBondInvariantsGenerator::use_bond_types,python:"use_bond_types",javascript:"useBondTypes",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganBondInvariantsGenerator)->bool,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganBondInvariantsGenerator.include_chirality",item:callable,owner:type_,rust:crate::MorganBondInvariantsGenerator::include_chirality,python:"include_chirality",javascript:"includeChirality",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganBondInvariantsGenerator)->bool,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganCallParams.new",item:callable,owner:type_,rust:crate::MorganCallParams::new,python:"new",javascript:"new",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:from_atoms,type:Option<Vec<u32>>,default:none},{name:ignore_atoms,type:Option<Vec<u32>>,default:none},{name:custom_atom_invariants,type:Option<Vec<u32>>,default:none},{name:custom_bond_invariants,type:Option<Vec<u32>>,default:none},{name:conformer_id,type:i32,default:"-1"}],output:crate::MorganCallParams,error:none,state:value_returning,operation:none,signature:fn(Option<Vec<u32>>,Option<Vec<u32>>,Option<Vec<u32>>,Option<Vec<u32>>,i32)->crate::MorganCallParams,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.new",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::new,python:"new",javascript:"new",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:params,type:Option<&'a crate::MorganParams>,default:none},{name:atom_invariants,type:Option<&'b crate::MorganAtomInvariantsGenerator>,default:none},{name:bond_invariants,type:Option<&'c crate::MorganBondInvariantsGenerator>,default:none}],output:crate::MorganFingerprintGenerator,error:crate::MorganReadError,state:value_returning,operation:none,signature:for<'a,'b,'c> fn(Option<&'a crate::MorganParams>,Option<&'b crate::MorganAtomInvariantsGenerator>,Option<&'c crate::MorganBondInvariantsGenerator>)->Result<crate::MorganFingerprintGenerator,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.from_json",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::from_json,python:"from_json",javascript:"fromJson",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:json,type:&'a str,default:required}],output:crate::MorganFingerprintGenerator,error:crate::MorganReadError,state:value_returning,operation:none,signature:for<'a> fn(&'a str)->Result<crate::MorganFingerprintGenerator,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.settings",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::settings,python:"settings",javascript:"settings",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:crate::MorganSettings,error:none,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator)->crate::MorganSettings,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.info_string",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::info_string,python:"info_string",javascript:"infoString",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:String,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator)->Result<String,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.to_json",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::to_json,python:"to_json",javascript:"toJson",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:String,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator)->Result<String,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.fingerprints",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::fingerprints,python:"fingerprints",javascript:"fingerprints",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:molecules,type:&'b [Option<&'b crate::Molecule>],default:required},{name:num_threads,type:i32,default:"1"}],output:Vec<Option<crate::Fingerprint>>,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::Fingerprint>>,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.morgan_fingerprint_with_generator",item:callable,owner:molecule,rust:crate::Molecule::morgan_fingerprint_with_generator,python:"morgan_fingerprint_with_generator",javascript:"morganFingerprintWithGenerator",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:generator,type:&'b crate::MorganFingerprintGenerator,default:required},{name:params,type:Option<&'c crate::MorganCallParams>,default:none},{name:output,type:Option<&'d mut crate::FingerprintAdditionalOutput>,default:none}],output:crate::Fingerprint,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::MorganFingerprintGenerator,Option<&'c crate::MorganCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::Fingerprint,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.counts",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::counts,python:"counts",javascript:"counts",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:molecules,type:&'b [Option<&'b crate::Molecule>],default:required},{name:num_threads,type:i32,default:"1"}],output:Vec<Option<crate::SparseCountFingerprint32>>,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseCountFingerprint32>>,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.morgan_count_fingerprint_with_generator",item:callable,owner:molecule,rust:crate::Molecule::morgan_count_fingerprint_with_generator,python:"morgan_count_fingerprint_with_generator",javascript:"morganCountFingerprintWithGenerator",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:generator,type:&'b crate::MorganFingerprintGenerator,default:required},{name:params,type:Option<&'c crate::MorganCallParams>,default:none},{name:output,type:Option<&'d mut crate::FingerprintAdditionalOutput>,default:none}],output:crate::SparseCountFingerprint32,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::MorganFingerprintGenerator,Option<&'c crate::MorganCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseCountFingerprint32,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.sparse_fingerprints",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::sparse_fingerprints,python:"sparse_fingerprints",javascript:"sparseFingerprints",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:molecules,type:&'b [Option<&'b crate::Molecule>],default:required},{name:num_threads,type:i32,default:"1"}],output:Vec<Option<crate::SparseBitFingerprint>>,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseBitFingerprint>>,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.morgan_sparse_fingerprint_with_generator",item:callable,owner:molecule,rust:crate::Molecule::morgan_sparse_fingerprint_with_generator,python:"morgan_sparse_fingerprint_with_generator",javascript:"morganSparseFingerprintWithGenerator",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:generator,type:&'b crate::MorganFingerprintGenerator,default:required},{name:params,type:Option<&'c crate::MorganCallParams>,default:none},{name:output,type:Option<&'d mut crate::FingerprintAdditionalOutput>,default:none}],output:crate::SparseBitFingerprint,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::MorganFingerprintGenerator,Option<&'c crate::MorganCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseBitFingerprint,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganFingerprintGenerator.sparse_counts",item:callable,owner:type_,rust:crate::MorganFingerprintGenerator::sparse_counts,python:"sparse_counts",javascript:"sparseCounts",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:molecules,type:&'b [Option<&'b crate::Molecule>],default:required},{name:num_threads,type:i32,default:"1"}],output:Vec<Option<crate::SparseCountFingerprint>>,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseCountFingerprint>>,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.morgan_sparse_count_fingerprint_with_generator",item:callable,owner:molecule,rust:crate::Molecule::morgan_sparse_count_fingerprint_with_generator,python:"morgan_sparse_count_fingerprint_with_generator",javascript:"morganSparseCountFingerprintWithGenerator",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[{name:generator,type:&'b crate::MorganFingerprintGenerator,default:required},{name:params,type:Option<&'c crate::MorganCallParams>,default:none},{name:output,type:Option<&'d mut crate::FingerprintAdditionalOutput>,default:none}],output:crate::SparseCountFingerprint,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::MorganFingerprintGenerator,Option<&'c crate::MorganCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseCountFingerprint,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.radius",item:callable,owner:type_,rust:crate::MorganSettings::radius,python:"radius",javascript:"radius",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:u32,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<u32,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_radius",item:callable,owner:type_,rust:crate::MorganSettings::set_radius,python:"set_radius",javascript:"setRadius",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:u32,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,u32)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.only_nonzero_invariants",item:callable,owner:type_,rust:crate::MorganSettings::only_nonzero_invariants,python:"only_nonzero_invariants",javascript:"onlyNonzeroInvariants",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<bool,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_only_nonzero_invariants",item:callable,owner:type_,rust:crate::MorganSettings::set_only_nonzero_invariants,python:"set_only_nonzero_invariants",javascript:"setOnlyNonzeroInvariants",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:bool,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,bool)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.include_redundant_environments",item:callable,owner:type_,rust:crate::MorganSettings::include_redundant_environments,python:"include_redundant_environments",javascript:"includeRedundantEnvironments",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<bool,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_include_redundant_environments",item:callable,owner:type_,rust:crate::MorganSettings::set_include_redundant_environments,python:"set_include_redundant_environments",javascript:"setIncludeRedundantEnvironments",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:bool,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,bool)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.include_chirality",item:callable,owner:type_,rust:crate::MorganSettings::include_chirality,python:"include_chirality",javascript:"includeChirality",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<bool,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_include_chirality",item:callable,owner:type_,rust:crate::MorganSettings::set_include_chirality,python:"set_include_chirality",javascript:"setIncludeChirality",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:bool,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,bool)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.count_simulation",item:callable,owner:type_,rust:crate::MorganSettings::count_simulation,python:"count_simulation",javascript:"countSimulation",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:bool,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<bool,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_count_simulation",item:callable,owner:type_,rust:crate::MorganSettings::set_count_simulation,python:"set_count_simulation",javascript:"setCountSimulation",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:bool,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,bool)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.fp_size",item:callable,owner:type_,rust:crate::MorganSettings::fp_size,python:"fp_size",javascript:"fpSize",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:u32,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<u32,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_fp_size",item:callable,owner:type_,rust:crate::MorganSettings::set_fp_size,python:"set_fp_size",javascript:"setFpSize",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:u32,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,u32)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.bits_per_feature",item:callable,owner:type_,rust:crate::MorganSettings::bits_per_feature,python:"bits_per_feature",javascript:"bitsPerFeature",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:u32,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<u32,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_bits_per_feature",item:callable,owner:type_,rust:crate::MorganSettings::set_bits_per_feature,python:"set_bits_per_feature",javascript:"setBitsPerFeature",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:u32,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,u32)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.count_bounds",item:callable,owner:type_,rust:crate::MorganSettings::count_bounds,python:"count_bounds",javascript:"countBounds",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:Vec<u32>,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<Vec<u32>,crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.set_count_bounds",item:callable,owner:type_,rust:crate::MorganSettings::set_count_bounds,python:"set_count_bounds",javascript:"setCountBounds",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:mutable,parameters:[{name:value,type:Vec<u32>,default:required}],output:(),error:crate::MorganReadError,state:in_place,operation:none,signature:for<'a,'b,'c,'d> fn(&'a mut crate::MorganSettings,Vec<u32>)->Result<(),crate::MorganReadError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"MorganSettings.params",item:callable,owner:type_,rust:crate::MorganSettings::params,python:"params",javascript:"params",feature:"cap-fingerprints",status:experimental,kind:instance,parameters:[],output:crate::MorganParams,error:crate::MorganReadError,state:read_only,operation:none,signature:for<'a,'b,'c,'d> fn(&'a crate::MorganSettings)->Result<crate::MorganParams,crate::MorganReadError>,},
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.FingerprintPreparationError", item: type, owner: type_,
            rust: crate::FingerprintPreparationError,
            python: "FingerprintPreparationError", javascript: "FingerprintPreparationError",
            feature: "cap-fingerprints", status: experimental, role: error,
        },
#[cfg(feature="cap-forcefields")]
{semantic_id:"types.UffEvaluationParams",item:type,owner:type_,rust:crate::UffEvaluationParams,python:"UffEvaluationParams",javascript:"UffEvaluationParams",feature:"cap-forcefields",status:experimental,role:parameter,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"types.UffEnergyGradient",item:type,owner:type_,rust:crate::UffEnergyGradient,python:"UffEnergyGradient",javascript:"UffEnergyGradient",feature:"cap-forcefields",status:experimental,role:result,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"Molecule.uff_energy_gradient",item:callable,owner:molecule,rust:crate::Molecule::uff_energy_gradient,python:"uff_energy_gradient",javascript:"uffEnergyGradient",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:crate::UffEnergyGradient,error:crate::OperationError,state:read_only,operation:none,signature:fn(&crate::Molecule)->Result<crate::UffEnergyGradient,crate::OperationError>,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"Molecule.uff_energy_gradient_with_params",item:callable,owner:molecule,rust:crate::Molecule::uff_energy_gradient_with_params,python:"uff_energy_gradient_with_params",javascript:"uffEnergyGradientWithParams",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[{name:params,type:&crate::UffEvaluationParams,default:required}],output:crate::UffEnergyGradient,error:crate::OperationError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::UffEvaluationParams)->Result<crate::UffEnergyGradient,crate::OperationError>,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"UffEnergyGradient.energy",item:callable,owner:type_,rust:crate::UffEnergyGradient::energy,python:"energy",javascript:"energy",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::UffEnergyGradient)->f64,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"UffEnergyGradient.gradient",item:callable,owner:type_,rust:crate::UffEnergyGradient::gradient,python:"gradient",javascript:"gradient",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a [f64],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::UffEnergyGradient)->&'a [f64],},

        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffOptimizationResult.molecule",item:callable,owner:type_,rust:crate::UffOptimizationResult::molecule,python:"molecule",javascript:"molecule",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a crate::Molecule,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::UffOptimizationResult)->&'a crate::Molecule,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffOptimizationResult.status_code",item:callable,owner:type_,rust:crate::UffOptimizationResult::status_code,python:"status_code",javascript:"statusCode",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:i32,error:none,state:read_only,operation:none,signature:fn(& crate::UffOptimizationResult)->i32,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffOptimizationResult.needs_more",item:callable,owner:type_,rust:crate::UffOptimizationResult::needs_more,python:"needs_more",javascript:"needsMore",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(& crate::UffOptimizationResult)->bool,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffOptimizationResult.energy",item:callable,owner:type_,rust:crate::UffOptimizationResult::energy,python:"energy",javascript:"energy",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(& crate::UffOptimizationResult)->f64,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerResult.conformer_id",item:callable,owner:type_,rust:crate::UffConformerResult::conformer_id,python:"conformer_id",javascript:"conformerId",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:usize,error:none,state:read_only,operation:none,signature:fn(& crate::UffConformerResult)->usize,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerResult.status_code",item:callable,owner:type_,rust:crate::UffConformerResult::status_code,python:"status_code",javascript:"statusCode",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:i32,error:none,state:read_only,operation:none,signature:fn(& crate::UffConformerResult)->i32,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerResult.needs_more",item:callable,owner:type_,rust:crate::UffConformerResult::needs_more,python:"needs_more",javascript:"needsMore",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(& crate::UffConformerResult)->bool,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerResult.energy",item:callable,owner:type_,rust:crate::UffConformerResult::energy,python:"energy",javascript:"energy",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(& crate::UffConformerResult)->f64,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerOptimizationResult.molecule",item:callable,owner:type_,rust:crate::UffConformerOptimizationResult::molecule,python:"molecule",javascript:"molecule",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a crate::Molecule,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::UffConformerOptimizationResult)->&'a crate::Molecule,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"UffConformerOptimizationResult.conformer_results",item:callable,owner:type_,rust:crate::UffConformerOptimizationResult::conformer_results,python:"conformer_results",javascript:"conformerResults",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a [crate::UffConformerResult],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::UffConformerOptimizationResult)->&'a [crate::UffConformerResult],},
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.AtomPairAtomCodeResult", item: type, owner: type_,
            rust: crate::AtomPairAtomCodeResult, python: "AtomPairAtomCodeResult", javascript: "AtomPairAtomCodeResult",
            feature: "cap-fingerprints", role: result,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.with_atom_pair_atom_code", item: callable, owner: molecule,
            rust: crate::Molecule::with_atom_pair_atom_code, python: "with_atom_pair_atom_code", javascript: "withAtomPairAtomCode",
            feature: "cap-fingerprints", kind: instance,
            parameters: [
                { name: atom_id, type: crate::AtomId, default: required },
                { name: branch_subtract, type: u32, default: "0" },
                { name: include_chirality, type: bool, default: "false" },
                { name: use_legacy_stereo_perception, type: bool, default: "true" },
            ],
            output: crate::AtomPairAtomCodeResult, error: crate::OperationError,
            state: value_returning, operation: "with_atom_pair_atom_code",
            signature: fn(&crate::Molecule, crate::AtomId, u32, bool, bool) -> Result<crate::AtomPairAtomCodeResult, crate::OperationError>,
        },

#[cfg(feature="cap-forcefields")]
{semantic_id:"MmffAtomProperties.atom_type",item:callable,owner:type_,rust:crate::MmffAtomProperties::atom_type,python:"atom_type",javascript:"atomType",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:u8,error:none,state:read_only,operation:none,signature:fn(&crate::MmffAtomProperties)->u8,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"MmffAtomProperties.formal_charge",item:callable,owner:type_,rust:crate::MmffAtomProperties::formal_charge,python:"formal_charge",javascript:"formalCharge",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::MmffAtomProperties)->f64,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"MmffAtomProperties.partial_charge",item:callable,owner:type_,rust:crate::MmffAtomProperties::partial_charge,python:"partial_charge",javascript:"partialCharge",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::MmffAtomProperties)->f64,},

#[cfg(feature="cap-forcefields")]
{semantic_id:"types.MmffEvaluationParams",item:type,owner:type_,rust:crate::MmffEvaluationParams,python:"MmffEvaluationParams",javascript:"MmffEvaluationParams",feature:"cap-forcefields",status:experimental,role:parameter,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"types.MmffEnergyGradient",item:type,owner:type_,rust:crate::MmffEnergyGradient,python:"MmffEnergyGradient",javascript:"MmffEnergyGradient",feature:"cap-forcefields",status:experimental,role:result,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"Molecule.mmff_energy_gradient",item:callable,owner:molecule,rust:crate::Molecule::mmff_energy_gradient,python:"mmff_energy_gradient",javascript:"mmffEnergyGradient",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:Option<crate::MmffEnergyGradient>,error:crate::OperationError,state:read_only,operation:none,signature:fn(&crate::Molecule)->Result<Option<crate::MmffEnergyGradient>,crate::OperationError>,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"Molecule.mmff_energy_gradient_with_params",item:callable,owner:molecule,rust:crate::Molecule::mmff_energy_gradient_with_params,python:"mmff_energy_gradient_with_params",javascript:"mmffEnergyGradientWithParams",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[{name:params,type:&crate::MmffEvaluationParams,default:required}],output:Option<crate::MmffEnergyGradient>,error:crate::OperationError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::MmffEvaluationParams)->Result<Option<crate::MmffEnergyGradient>,crate::OperationError>,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"MmffEnergyGradient.energy",item:callable,owner:type_,rust:crate::MmffEnergyGradient::energy,python:"energy",javascript:"energy",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::MmffEnergyGradient)->f64,},
#[cfg(feature="cap-forcefields")]
{semantic_id:"MmffEnergyGradient.gradient",item:callable,owner:type_,rust:crate::MmffEnergyGradient::gradient,python:"gradient",javascript:"gradient",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a [f64],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MmffEnergyGradient)->&'a [f64],},

        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeResult.molecule",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeResult::molecule,python:"molecule",javascript:"molecule",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a crate::Molecule,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MmffOptimizeMoleculeResult)->&'a crate::Molecule,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeResult.needs_more",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeResult::needs_more,python:"needs_more",javascript:"needsMore",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::MmffOptimizeMoleculeResult)->bool,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeResult.status_code",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeResult::status_code,python:"status_code",javascript:"statusCode",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:i32,error:none,state:read_only,operation:none,signature:fn(&crate::MmffOptimizeMoleculeResult)->i32,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeConfsResult.molecule",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeConfsResult::molecule,python:"molecule",javascript:"molecule",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a crate::Molecule,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MmffOptimizeMoleculeConfsResult)->&'a crate::Molecule,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeConfsResult.conformer_results",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeConfsResult::conformer_results,python:"conformer_results",javascript:"conformerResults",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a [crate::MmffOptimizeMoleculeConfResult],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MmffOptimizeMoleculeConfsResult)->&'a [crate::MmffOptimizeMoleculeConfResult],},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeConfResult.needs_more",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeConfResult::needs_more,python:"needs_more",javascript:"needsMore",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::MmffOptimizeMoleculeConfResult)->bool,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeConfResult.status_code",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeConfResult::status_code,python:"status_code",javascript:"statusCode",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:i32,error:none,state:read_only,operation:none,signature:fn(&crate::MmffOptimizeMoleculeConfResult)->i32,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffOptimizeMoleculeConfResult.energy",item:callable,owner:type_,rust:crate::MmffOptimizeMoleculeConfResult::energy,python:"energy",javascript:"energy",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::MmffOptimizeMoleculeConfResult)->f64,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.is_valid",item:callable,owner:type_,rust:crate::MmffProperties::is_valid,python:"is_valid",javascript:"isValid",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::MmffProperties)->bool,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.variant",item:callable,owner:type_,rust:crate::MmffProperties::variant,python:"variant",javascript:"variant",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:crate::MmffVariant,error:none,state:read_only,operation:none,signature:fn(&crate::MmffProperties)->crate::MmffVariant,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.atoms",item:callable,owner:type_,rust:crate::MmffProperties::atoms,python:"atoms",javascript:"atoms",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[],output:&'a [crate::MmffAtomProperties],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MmffProperties)->&'a [crate::MmffAtomProperties],},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.atom_type",item:callable,owner:type_,rust:crate::MmffProperties::atom_type,python:"atom_type",javascript:"atomType",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[{ name: atom_index, type: usize, default: required }],output:u8,error:crate::MmffMolPropertiesError,state:read_only,operation:none,signature:fn(&crate::MmffProperties,usize)->Result<u8,crate::MmffMolPropertiesError>,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.formal_charge",item:callable,owner:type_,rust:crate::MmffProperties::formal_charge,python:"formal_charge",javascript:"formalCharge",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[{ name: atom_index, type: usize, default: required }],output:f64,error:crate::MmffMolPropertiesError,state:read_only,operation:none,signature:fn(&crate::MmffProperties,usize)->Result<f64,crate::MmffMolPropertiesError>,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"MmffProperties.partial_charge",item:callable,owner:type_,rust:crate::MmffProperties::partial_charge,python:"partial_charge",javascript:"partialCharge",feature:"cap-forcefields",status:experimental,kind:instance,parameters:[{ name: atom_index, type: usize, default: required }],output:f64,error:crate::MmffMolPropertiesError,state:read_only,operation:none,signature:fn(&crate::MmffProperties,usize)->Result<f64,crate::MmffMolPropertiesError>,},

        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffOptimizationParams", item: type, owner: type_, rust: crate::MmffOptimizationParams, python: "MmffOptimizationParams", javascript: "MmffOptimizationParams", feature: "cap-forcefields", status: experimental, role: parameter, },
        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffConformerOptimizationParams", item: type, owner: type_, rust: crate::MmffConformerOptimizationParams, python: "MmffConformerOptimizationParams", javascript: "MmffConformerOptimizationParams", feature: "cap-forcefields", status: experimental, role: parameter, },
        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffOptimizeMoleculeResult", item: type, owner: type_, rust: crate::MmffOptimizeMoleculeResult, python: "MmffOptimizeMoleculeResult", javascript: "MmffOptimizeMoleculeResult", feature: "cap-forcefields", status: experimental, role: result, },
        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffOptimizeMoleculeConfResult", item: type, owner: type_, rust: crate::MmffOptimizeMoleculeConfResult, python: "MmffOptimizeMoleculeConfResult", javascript: "MmffOptimizeMoleculeConfResult", feature: "cap-forcefields", status: experimental, role: result, },
        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffOptimizeMoleculeConfsResult", item: type, owner: type_, rust: crate::MmffOptimizeMoleculeConfsResult, python: "MmffOptimizeMoleculeConfsResult", javascript: "MmffOptimizeMoleculeConfsResult", feature: "cap-forcefields", status: experimental, role: result, },
        #[cfg(feature="cap-forcefields")]
        { semantic_id: "types.MmffOptimizationError", item: type, owner: type_, rust: crate::MmffOptimizationError, python: "MmffOptimizationError", javascript: "MmffOptimizationError", feature: "cap-forcefields", status: experimental, role: error, },
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"Molecule.with_mmff_optimized",item:callable,owner:molecule,rust:crate::Molecule::with_mmff_optimized,python:"with_mmff_optimized",javascript:"withMmffOptimized",feature:"cap-forcefields",kind:instance,parameters:[],output:crate::MmffOptimizeMoleculeResult,error:crate::OperationError,state:value_returning,operation:"with_mmff_optimized",signature:fn(&crate::Molecule) -> Result<crate::MmffOptimizeMoleculeResult, crate::OperationError>,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"Molecule.with_mmff_optimized_with_params",item:callable,owner:molecule,rust:crate::Molecule::with_mmff_optimized_with_params,python:"with_mmff_optimized_with_params",javascript:"withMmffOptimizedWithParams",feature:"cap-forcefields",kind:instance,parameters:[{ name: params, type: &crate::MmffOptimizationParams, default: required }],output:crate::MmffOptimizeMoleculeResult,error:crate::OperationError,state:value_returning,operation:"with_mmff_optimized_with_params",signature:fn(&crate::Molecule, &crate::MmffOptimizationParams) -> Result<crate::MmffOptimizeMoleculeResult, crate::OperationError>,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"Molecule.with_mmff_optimized_confs",item:callable,owner:molecule,rust:crate::Molecule::with_mmff_optimized_confs,python:"with_mmff_optimized_confs",javascript:"withMmffOptimizedConfs",feature:"cap-forcefields",kind:instance,parameters:[],output:crate::MmffOptimizeMoleculeConfsResult,error:crate::OperationError,state:value_returning,operation:"with_mmff_optimized_confs",signature:fn(&crate::Molecule) -> Result<crate::MmffOptimizeMoleculeConfsResult, crate::OperationError>,},
        #[cfg(feature="cap-forcefields")]
        {semantic_id:"Molecule.with_mmff_optimized_confs_with_params",item:callable,owner:molecule,rust:crate::Molecule::with_mmff_optimized_confs_with_params,python:"with_mmff_optimized_confs_with_params",javascript:"withMmffOptimizedConfsWithParams",feature:"cap-forcefields",kind:instance,parameters:[{ name: params, type: &crate::MmffConformerOptimizationParams, default: required }],output:crate::MmffOptimizeMoleculeConfsResult,error:crate::OperationError,state:value_returning,operation:"with_mmff_optimized_confs_with_params",signature:fn(&crate::Molecule, &crate::MmffConformerOptimizationParams) -> Result<crate::MmffOptimizeMoleculeConfsResult, crate::OperationError>,},

        #[cfg(feature="cap-tautomer")]
        { semantic_id:"default_tautomer_score_terms", item:callable, owner:module, rust:crate::default_tautomer_score_terms, python:"default_tautomer_score_terms", javascript:"defaultTautomerScoreTerms", feature:"cap-tautomer", status:experimental, kind:module,
          parameters:[], output:&'static [crate::TautomerScoreTerm], error:none, state:value_returning, operation:none,
          signature: fn() -> &'static [crate::TautomerScoreTerm], },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.max_tautomers", item:callable, owner:type_, rust:crate::TautomerParams::max_tautomers, python:"max_tautomers", javascript:"maxTautomers", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:u32, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> u32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_max_tautomers", item:callable, owner:type_, rust:crate::TautomerParams::set_max_tautomers, python:"set_max_tautomers", javascript:"setMaxTautomers", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:u32,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, u32) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_max_tautomers", item:callable, owner:type_, rust:crate::TautomerParams::with_max_tautomers, python:"with_max_tautomers", javascript:"withMaxTautomers", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:u32,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, u32) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.max_transforms", item:callable, owner:type_, rust:crate::TautomerParams::max_transforms, python:"max_transforms", javascript:"maxTransforms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:u32, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> u32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_max_transforms", item:callable, owner:type_, rust:crate::TautomerParams::set_max_transforms, python:"set_max_transforms", javascript:"setMaxTransforms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:u32,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, u32) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_max_transforms", item:callable, owner:type_, rust:crate::TautomerParams::with_max_transforms, python:"with_max_transforms", javascript:"withMaxTransforms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:u32,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, u32) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.remove_sp3_stereo", item:callable, owner:type_, rust:crate::TautomerParams::remove_sp3_stereo, python:"remove_sp3_stereo", javascript:"removeSp3Stereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_remove_sp3_stereo", item:callable, owner:type_, rust:crate::TautomerParams::set_remove_sp3_stereo, python:"set_remove_sp3_stereo", javascript:"setRemoveSp3Stereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:bool,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, bool) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_remove_sp3_stereo", item:callable, owner:type_, rust:crate::TautomerParams::with_remove_sp3_stereo, python:"with_remove_sp3_stereo", javascript:"withRemoveSp3Stereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:bool,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, bool) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.remove_bond_stereo", item:callable, owner:type_, rust:crate::TautomerParams::remove_bond_stereo, python:"remove_bond_stereo", javascript:"removeBondStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_remove_bond_stereo", item:callable, owner:type_, rust:crate::TautomerParams::set_remove_bond_stereo, python:"set_remove_bond_stereo", javascript:"setRemoveBondStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:bool,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, bool) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_remove_bond_stereo", item:callable, owner:type_, rust:crate::TautomerParams::with_remove_bond_stereo, python:"with_remove_bond_stereo", javascript:"withRemoveBondStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:bool,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, bool) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.remove_isotopic_hydrogens", item:callable, owner:type_, rust:crate::TautomerParams::remove_isotopic_hydrogens, python:"remove_isotopic_hydrogens", javascript:"removeIsotopicHydrogens", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_remove_isotopic_hydrogens", item:callable, owner:type_, rust:crate::TautomerParams::set_remove_isotopic_hydrogens, python:"set_remove_isotopic_hydrogens", javascript:"setRemoveIsotopicHydrogens", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:bool,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, bool) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_remove_isotopic_hydrogens", item:callable, owner:type_, rust:crate::TautomerParams::with_remove_isotopic_hydrogens, python:"with_remove_isotopic_hydrogens", javascript:"withRemoveIsotopicHydrogens", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:bool,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, bool) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.reassign_stereo", item:callable, owner:type_, rust:crate::TautomerParams::reassign_stereo, python:"reassign_stereo", javascript:"reassignStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_reassign_stereo", item:callable, owner:type_, rust:crate::TautomerParams::set_reassign_stereo, python:"set_reassign_stereo", javascript:"setReassignStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:value,type:bool,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, bool) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.with_reassign_stereo", item:callable, owner:type_, rust:crate::TautomerParams::with_reassign_stereo, python:"with_reassign_stereo", javascript:"withReassignStereo", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[{name:value,type:bool,default:required}], output:crate::TautomerParams, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerParams, bool) -> crate::TautomerParams, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.v1", item:callable, owner:type_, rust:crate::TautomerParams::v1, python:"v1", javascript:"v1", feature:"cap-tautomer", status:experimental, kind:static_,
          parameters:[], output:crate::TautomerParams, error:crate::TautomerCatalogError, state:value_returning, operation:none,
          signature: fn() -> Result<crate::TautomerParams,crate::TautomerCatalogError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.from_transform_data", item:callable, owner:type_, rust:crate::TautomerParams::from_transform_data, python:"from_transform_data", javascript:"fromTransformData", feature:"cap-tautomer", status:experimental, kind:static_,
          parameters:[{name:data,type:&'a [(&'b str,&'b str,&'b str,&'b str)],default:required}], output:crate::TautomerParams, error:crate::TautomerCatalogError, state:value_returning, operation:none,
          signature:for<'a,'b:'a> fn(&'a [(&'b str,&'b str,&'b str,&'b str)]) -> Result<crate::TautomerParams,crate::TautomerCatalogError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.from_transform_file", item:callable, owner:type_, rust:crate::TautomerParams::from_transform_file, python:"from_transform_file", javascript:"fromTransformFile", feature:"cap-tautomer", status:experimental, kind:static_,
          parameters:[{name:path,type:&'a str,default:required}], output:crate::TautomerParams, error:crate::TautomerCatalogError, state:value_returning, operation:none,
          signature:for<'a> fn(&'a str) -> Result<crate::TautomerParams,crate::TautomerCatalogError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.transform_count", item:callable, owner:type_, rust:crate::TautomerParams::transform_count, python:"transform_count", javascript:"transformCount", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:usize, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> usize, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.callback", item:callable, owner:type_, rust:crate::TautomerParams::callback, python:"callback", javascript:"callback", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:Option<&'a std::sync::Arc<dyn crate::TautomerEnumerationCallback>>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> Option<&'a std::sync::Arc<dyn crate::TautomerEnumerationCallback>>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_callback", item:callable, owner:type_, rust:crate::TautomerParams::set_callback, python:"set_callback", javascript:"setCallback", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:callback,type:Option<std::sync::Arc<dyn crate::TautomerEnumerationCallback>>,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, Option<std::sync::Arc<dyn crate::TautomerEnumerationCallback>>) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.scorer", item:callable, owner:type_, rust:crate::TautomerParams::scorer, python:"scorer", javascript:"scorer", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:Option<&'a std::sync::Arc<dyn crate::TautomerScorer>>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerParams) -> Option<&'a std::sync::Arc<dyn crate::TautomerScorer>>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerParams.set_scorer", item:callable, owner:type_, rust:crate::TautomerParams::set_scorer, python:"set_scorer", javascript:"setScorer", feature:"cap-tautomer", status:experimental, kind:instance, receiver:mutable,
          parameters:[{name:scorer,type:Option<std::sync::Arc<dyn crate::TautomerScorer>>,default:required}], output:(), error:none, state:in_place, operation:none,
          signature:for<'a> fn(&'a mut crate::TautomerParams, Option<std::sync::Arc<dyn crate::TautomerScorer>>) -> (), },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.len", item:callable, owner:type_, rust:crate::TautomerEnumeration::len, python:"len", javascript:"len", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:usize, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> usize, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.is_empty", item:callable, owner:type_, rust:crate::TautomerEnumeration::is_empty, python:"is_empty", javascript:"isEmpty", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.status", item:callable, owner:type_, rust:crate::TautomerEnumeration::status, python:"status", javascript:"status", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:crate::TautomerEnumerationStatus, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> crate::TautomerEnumerationStatus, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.modified_atoms", item:callable, owner:type_, rust:crate::TautomerEnumeration::modified_atoms, python:"modified_atoms", javascript:"modifiedAtoms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a std::collections::BTreeSet<crate::AtomId>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> &'a std::collections::BTreeSet<crate::AtomId>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.modified_bonds", item:callable, owner:type_, rust:crate::TautomerEnumeration::modified_bonds, python:"modified_bonds", javascript:"modifiedBonds", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a std::collections::BTreeSet<crate::BondId>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> &'a std::collections::BTreeSet<crate::BondId>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.canonical_smiles", item:callable, owner:type_, rust:crate::TautomerEnumeration::canonical_smiles, python:"canonical_smiles", javascript:"canonicalSmiles", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:Vec<&'a str>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> Vec<&'a str>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.get", item:callable, owner:type_, rust:crate::TautomerEnumeration::get, python:"get", javascript:"get", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[{name:index,type:usize,default:required}], output:Option<&'a crate::Molecule>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration, usize) -> Option<&'a crate::Molecule>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.iter", item:callable, owner:type_, rust:crate::TautomerEnumeration::iter, python:"iter", javascript:"iter", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:std::iter::Map<std::slice::Iter<'a,(String,crate::Molecule)>,for<'b> fn(&'b (String,crate::Molecule))-> &'b crate::Molecule>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> std::iter::Map<std::slice::Iter<'a,(String,crate::Molecule)>,for<'b> fn(&'b (String,crate::Molecule))-> &'b crate::Molecule>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.entries", item:callable, owner:type_, rust:crate::TautomerEnumeration::entries, python:"entries", javascript:"entries", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:std::iter::Map<std::slice::Iter<'a,(String,crate::Molecule)>,for<'b> fn(&'b (String,crate::Molecule))-> (&'b str,&'b crate::Molecule)>, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerEnumeration) -> std::iter::Map<std::slice::Iter<'a,(String,crate::Molecule)>,for<'b> fn(&'b (String,crate::Molecule))-> (&'b str,&'b crate::Molecule)>, },
        #[cfg(feature="cap-tautomer")]
        {semantic_id:"tautomer.canonical_tautomer_from_molecules",item:callable,owner:module,rust:crate::canonical_tautomer_from_molecules,python:"canonical_tautomer_from_molecules",javascript:"canonicalTautomerFromMolecules",feature:"cap-tautomer",status:experimental,kind:module,
          parameters:[{name:molecules,type:&[crate::Molecule],default:required}],output:crate::Molecule,error:crate::OperationError,state:value_returning,operation:none,
          signature:for<'a> fn(&'a [crate::Molecule]) -> Result<crate::Molecule,crate::OperationError>,},
        #[cfg(feature="cap-tautomer")]
        {semantic_id:"tautomer.canonical_tautomer_from_molecules_with_params",item:callable,owner:module,rust:crate::canonical_tautomer_from_molecules_with_params,python:"canonical_tautomer_from_molecules_with_params",javascript:"canonicalTautomerFromMoleculesWithParams",feature:"cap-tautomer",status:experimental,kind:module,
          parameters:[{name:molecules,type:&[crate::Molecule],default:required},{name:params,type:&crate::TautomerParams,default:required}],output:crate::Molecule,error:crate::OperationError,state:value_returning,operation:none,
          signature:for<'a,'b> fn(&'a [crate::Molecule], &'b crate::TautomerParams) -> Result<crate::Molecule,crate::OperationError>,},
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.to_owned", item:callable, owner:type_, rust:crate::TautomerMoleculeView::to_owned, python:"to_owned", javascript:"toOwned", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:crate::TautomerMoleculeView<'static>, error:none, state:value_returning, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> crate::TautomerMoleculeView<'static>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.num_atoms", item:callable, owner:type_, rust:crate::TautomerMoleculeView::num_atoms, python:"num_atoms", javascript:"numAtoms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:usize, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> usize, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.num_bonds", item:callable, owner:type_, rust:crate::TautomerMoleculeView::num_bonds, python:"num_bonds", javascript:"numBonds", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:usize, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> usize, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.atom_metadata", item:callable, owner:type_, rust:crate::TautomerMoleculeView::atom_metadata, python:"atom_metadata", javascript:"atomMetadata", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:Vec<crate::AtomMetadata>, error:crate::ValenceError, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> Result<Vec<crate::AtomMetadata>,crate::ValenceError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.atom_degree", item:callable, owner:type_, rust:crate::TautomerMoleculeView::atom_degree, python:"atom_degree", javascript:"atomDegree", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[{name:id,type:crate::AtomId,default:required}], output:Option<usize>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>,crate::AtomId) -> Option<usize>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.atoms", item:callable, owner:type_, rust:crate::TautomerMoleculeView::atoms, python:"atoms", javascript:"atoms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a [crate::Atom], error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> &'a [crate::Atom], },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.bonds", item:callable, owner:type_, rust:crate::TautomerMoleculeView::bonds, python:"bonds", javascript:"bonds", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a [crate::Bond], error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> &'a [crate::Bond], },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.properties", item:callable, owner:type_, rust:crate::TautomerMoleculeView::properties, python:"properties", javascript:"properties", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a crate::MoleculeProperties, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> &'a crate::MoleculeProperties, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.atom", item:callable, owner:type_, rust:crate::TautomerMoleculeView::atom, python:"atom", javascript:"atom", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[{name:id,type:crate::AtomId,default:required}], output:Option<&'a crate::Atom>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>, crate::AtomId) -> Option<&'a crate::Atom>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.bond", item:callable, owner:type_, rust:crate::TautomerMoleculeView::bond, python:"bond", javascript:"bond", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[{name:id,type:crate::BondId,default:required}], output:Option<&'a crate::Bond>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>, crate::BondId) -> Option<&'a crate::Bond>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.to_smiles", item:callable, owner:type_, rust:crate::TautomerMoleculeView::to_smiles, python:"to_smiles", javascript:"toSmiles", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:String, error:crate::TautomerRunError, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> Result<String,crate::TautomerRunError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerMoleculeView.tautomer_score", item:callable, owner:type_, rust:crate::TautomerMoleculeView::tautomer_score, python:"tautomer_score", javascript:"tautomerScore", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:crate::TautomerScore, error:crate::TautomerRunError, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerMoleculeView<'b>) -> Result<crate::TautomerScore,crate::TautomerRunError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.to_owned", item:callable, owner:type_, rust:crate::TautomerProgress::to_owned, python:"to_owned", javascript:"toOwned", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:crate::TautomerProgress<'static>, error:none, state:value_returning, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> crate::TautomerProgress<'static>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.len", item:callable, owner:type_, rust:crate::TautomerProgress::len, python:"len", javascript:"len", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:usize, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> usize, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.is_empty", item:callable, owner:type_, rust:crate::TautomerProgress::is_empty, python:"is_empty", javascript:"isEmpty", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:bool, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> bool, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.status", item:callable, owner:type_, rust:crate::TautomerProgress::status, python:"status", javascript:"status", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:crate::TautomerEnumerationStatus, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> crate::TautomerEnumerationStatus, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.num_transforms", item:callable, owner:type_, rust:crate::TautomerProgress::num_transforms, python:"num_transforms", javascript:"numTransforms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:u32, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> u32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.modified_atoms", item:callable, owner:type_, rust:crate::TautomerProgress::modified_atoms, python:"modified_atoms", javascript:"modifiedAtoms", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a std::collections::BTreeSet<crate::AtomId>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> &'a std::collections::BTreeSet<crate::AtomId>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.modified_bonds", item:callable, owner:type_, rust:crate::TautomerProgress::modified_bonds, python:"modified_bonds", javascript:"modifiedBonds", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a std::collections::BTreeSet<crate::BondId>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> &'a std::collections::BTreeSet<crate::BondId>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerProgress.entries", item:callable, owner:type_, rust:crate::TautomerProgress::entries, python:"entries", javascript:"entries", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:Box<dyn ExactSizeIterator<Item=(&'a str,crate::TautomerMoleculeView<'a>)>+'a>, error:none, state:read_only, operation:none,
          signature:for<'a,'b:'a> fn(&'a crate::TautomerProgress<'b>) -> Box<dyn ExactSizeIterator<Item=(&'a str,crate::TautomerMoleculeView<'a>)>+'a>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScoreTerm.new", item:callable, owner:type_, rust:crate::TautomerScoreTerm::new, python:"new", javascript:"new", feature:"cap-tautomer", status:experimental, kind:static_,
          parameters:[{name:name,type:String,default:required}, {name:smarts,type:String,default:required}, {name:score,type:i32,default:required}], output:crate::TautomerScoreTerm, error:none, state:value_returning, operation:none,
          signature: fn(String, String, i32) -> crate::TautomerScoreTerm, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScoreTerm.name", item:callable, owner:type_, rust:crate::TautomerScoreTerm::name, python:"name", javascript:"name", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a str, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerScoreTerm) -> &'a str, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScoreTerm.smarts", item:callable, owner:type_, rust:crate::TautomerScoreTerm::smarts, python:"smarts", javascript:"smarts", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:&'a str, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerScoreTerm) -> &'a str, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScoreTerm.score", item:callable, owner:type_, rust:crate::TautomerScoreTerm::score, python:"score", javascript:"score", feature:"cap-tautomer", status:experimental, kind:instance, receiver:shared,
          parameters:[], output:i32, error:none, state:read_only, operation:none,
          signature:for<'a> fn(&'a crate::TautomerScoreTerm) -> i32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScore.ring", item:callable, owner:type_, rust:crate::TautomerScore::ring, python:"ring", javascript:"ring", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[], output:i32, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerScore) -> i32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScore.substructure", item:callable, owner:type_, rust:crate::TautomerScore::substructure, python:"substructure", javascript:"substructure", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[], output:i32, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerScore) -> i32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScore.hetero_hydrogen", item:callable, owner:type_, rust:crate::TautomerScore::hetero_hydrogen, python:"hetero_hydrogen", javascript:"heteroHydrogen", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[], output:i32, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerScore) -> i32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerScore.total", item:callable, owner:type_, rust:crate::TautomerScore::total, python:"total", javascript:"total", feature:"cap-tautomer", status:experimental, kind:instance, receiver:owned,
          parameters:[], output:i32, error:none, state:value_returning, operation:none,
          signature: fn(crate::TautomerScore) -> i32, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.canonical_tautomer", item:callable, owner:type_, rust:crate::TautomerEnumeration::canonical_tautomer, python:"canonical_tautomer", javascript:"canonicalTautomer", feature:"cap-tautomer", status:experimental, kind:instance,
           parameters:[], output:crate::Molecule, error:crate::OperationError, state:value_returning, operation:none,
           signature:for<'a> fn(&'a crate::TautomerEnumeration)->Result<crate::Molecule,crate::OperationError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"TautomerEnumeration.canonical_tautomer_with_params", item:callable, owner:type_, rust:crate::TautomerEnumeration::canonical_tautomer_with_params, python:"canonical_tautomer_with_params", javascript:"canonicalTautomerWithParams", feature:"cap-tautomer", status:experimental, kind:instance,
           parameters:[{name:params,type:&crate::TautomerParams,default:required}], output:crate::Molecule, error:crate::OperationError, state:value_returning, operation:none,
           signature:for<'a,'b> fn(&'a crate::TautomerEnumeration,&'b crate::TautomerParams)->Result<crate::Molecule,crate::OperationError>, },

        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerParams", item:type, owner:type_, rust:crate::TautomerParams, python:"TautomerParams", javascript:"TautomerParams", feature:"cap-tautomer", status:experimental, role:parameter, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerScoreParams", item:type, owner:type_, rust:crate::TautomerScoreParams, python:"TautomerScoreParams", javascript:"TautomerScoreParams", feature:"cap-tautomer", status:experimental, role:parameter, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerScoreTerm", item:type, owner:type_, rust:crate::TautomerScoreTerm, python:"TautomerScoreTerm", javascript:"TautomerScoreTerm", feature:"cap-tautomer", status:experimental, role:value, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerScore", item:type, owner:type_, rust:crate::TautomerScore, python:"TautomerScore", javascript:"TautomerScore", feature:"cap-tautomer", status:experimental, role:result, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerEnumeration", item:type, owner:type_, rust:crate::TautomerEnumeration, python:"TautomerEnumeration", javascript:"TautomerEnumeration", feature:"cap-tautomer", status:experimental, role:result, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerEnumerationStatus", item:type, owner:type_, rust:crate::TautomerEnumerationStatus, python:"TautomerEnumerationStatus", javascript:"TautomerEnumerationStatus", feature:"cap-tautomer", status:experimental, role:value, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerRunError", item:type, owner:type_, rust:crate::TautomerRunError, python:"TautomerRunError", javascript:"TautomerRunError", feature:"cap-tautomer", status:experimental, role:error, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerCatalogError", item:type, owner:type_, rust:crate::TautomerCatalogError, python:"TautomerCatalogError", javascript:"TautomerCatalogError", feature:"cap-tautomer", status:experimental, role:error, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerMoleculeView", item:type, owner:type_, rust:crate::TautomerMoleculeView, python:"TautomerMoleculeView", javascript:"TautomerMoleculeView", feature:"cap-tautomer", status:experimental, role:value, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"types.TautomerProgress", item:type, owner:type_, rust:crate::TautomerProgress, python:"TautomerProgress", javascript:"TautomerProgress", feature:"cap-tautomer", status:experimental, role:result, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.enumerate_tautomers", item:callable, owner:molecule, rust:crate::Molecule::enumerate_tautomers, python:"enumerate_tautomers", javascript:"enumerateTautomers", feature:"cap-tautomer", kind:instance,
           parameters:[], output:crate::TautomerEnumeration, error:crate::OperationError, state:value_returning, operation:"enumerate_tautomers",
           signature:for<'a> fn(&'a crate::Molecule)->Result<crate::TautomerEnumeration,crate::OperationError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.enumerate_tautomers_with_params", item:callable, owner:molecule, rust:crate::Molecule::enumerate_tautomers_with_params, python:"enumerate_tautomers_with_params", javascript:"enumerateTautomersWithParams", feature:"cap-tautomer", kind:instance,
           parameters:[{name:params,type:&crate::TautomerParams,default:required}], output:crate::TautomerEnumeration, error:crate::OperationError, state:value_returning, operation:"enumerate_tautomers_with_params",
           signature:for<'a, 'b> fn(&'a crate::Molecule, &'b crate::TautomerParams)->Result<crate::TautomerEnumeration,crate::OperationError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.canonical_tautomer", item:callable, owner:molecule, rust:crate::Molecule::canonical_tautomer, python:"canonical_tautomer", javascript:"canonicalTautomer", feature:"cap-tautomer", kind:instance,
           parameters:[], output:crate::Molecule, error:crate::OperationError, state:value_returning, operation:"canonical_tautomer",
           signature:for<'a> fn(&'a crate::Molecule)->Result<crate::Molecule,crate::OperationError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.canonical_tautomer_with_params", item:callable, owner:molecule, rust:crate::Molecule::canonical_tautomer_with_params, python:"canonical_tautomer_with_params", javascript:"canonicalTautomerWithParams", feature:"cap-tautomer", kind:instance,
           parameters:[{name:params,type:&crate::TautomerParams,default:required}], output:crate::Molecule, error:crate::OperationError, state:value_returning, operation:"canonical_tautomer_with_params",
           signature:for<'a, 'b> fn(&'a crate::Molecule, &'b crate::TautomerParams)->Result<crate::Molecule,crate::OperationError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.tautomer_score", item:callable, owner:molecule, rust:crate::Molecule::tautomer_score, python:"tautomer_score", javascript:"tautomerScore", feature:"cap-tautomer", status:experimental, kind:instance,
           parameters:[], output:crate::TautomerScore, error:crate::TautomerRunError, state:read_only, operation:none,
           signature:for<'a> fn(&'a crate::Molecule)->Result<crate::TautomerScore,crate::TautomerRunError>, },
        #[cfg(feature="cap-tautomer")]
        { semantic_id:"Molecule.tautomer_score_with_params", item:callable, owner:molecule, rust:crate::Molecule::tautomer_score_with_params, python:"tautomer_score_with_params", javascript:"tautomerScoreWithParams", feature:"cap-tautomer", status:experimental, kind:instance,
           parameters:[{name:params,type:&crate::TautomerScoreParams,default:required}], output:crate::TautomerScore, error:crate::TautomerRunError, state:read_only, operation:none,
           signature:for<'a, 'b> fn(&'a crate::Molecule, &'b crate::TautomerScoreParams)->Result<crate::TautomerScore,crate::TautomerRunError>, },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SmartsParseParams", item: type, owner: type_,
            rust: crate::SmartsParseParams, python: "SmartsParseParams", javascript: "SmartsParseParams",
            feature: "cap-search", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SmartsParseError", item: type, owner: type_,
            rust: crate::SmartsParseError, python: "SmartsParseError", javascript: "SmartsParseError",
            feature: "cap-search", status: experimental, role: error,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SubstructMatchParams", item: type, owner: type_,
            rust: crate::SubstructMatchParams, python: "SubstructMatchParams", javascript: "SubstructMatchParams",
            feature: "cap-search", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SubstructMatchError", item: type, owner: type_,
            rust: crate::SubstructMatchError, python: "SubstructMatchError", javascript: "SubstructMatchError",
            feature: "cap-search", status: experimental, role: error,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.MatchResult", item: type, owner: type_,
            rust: crate::MatchResult, python: "MatchResult", javascript: "MatchResult",
            feature: "cap-search", status: experimental, role: result,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.CompiledQuery", item: type, owner: type_,
            rust: crate::CompiledQuery, python: "CompiledQuery", javascript: "CompiledQuery",
            feature: "cap-search", status: experimental, role: value,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.QueryCompileError", item: type, owner: type_,
            rust: crate::QueryCompileError, python: "QueryCompileError", javascript: "QueryCompileError",
            feature: "cap-search", status: experimental, role: error,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.MatchError", item: type, owner: type_,
            rust: crate::MatchError, python: "MatchError", javascript: "MatchError",
            feature: "cap-search", status: experimental, role: error,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SmartsWriteParams", item: type, owner: type_,
            rust: crate::SmartsWriteParams, python: "SmartsWriteParams", javascript: "SmartsWriteParams",
            feature: "cap-search", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "types.SmartsWriteError", item: type, owner: type_,
            rust: crate::SmartsWriteError, python: "SmartsWriteError", javascript: "SmartsWriteError",
            feature: "cap-search", status: experimental, role: error,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "search.parse_smarts", item: callable, owner: module,
            rust: crate::parse_smarts, python: "parse_smarts", javascript: "parseSmarts",
            feature: "cap-search", status: experimental, kind: module,
            parameters: [{ name: text, type: &str, default: required }], output: crate::QueryGraph, error: crate::SmartsParseError,
            state: value_returning, operation: none,
            signature: for<'a> fn(&'a str) -> Result<crate::QueryGraph, crate::SmartsParseError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "search.parse_smarts_with_params", item: callable, owner: module,
            rust: crate::parse_smarts_with_params, python: "parse_smarts_with_params", javascript: "parseSmartsWithParams",
            feature: "cap-search", status: experimental, kind: module,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::SmartsParseParams, default: required }], output: crate::QueryGraph, error: crate::SmartsParseError,
            state: value_returning, operation: none,
            signature: for<'a, 'b> fn(&'a str, &'b crate::SmartsParseParams) -> Result<crate::QueryGraph, crate::SmartsParseError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "QueryGraph.from_smarts", item: callable, owner: type_,
            rust: crate::search::from_smarts, python: "from_smarts", javascript: "fromSmarts",
            feature: "cap-search", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }], output: crate::QueryGraph, error: crate::SmartsParseError,
            state: value_returning, operation: none,
            signature: for<'a> fn(&'a str) -> Result<crate::QueryGraph, crate::SmartsParseError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "QueryGraph.from_smarts_with_params", item: callable, owner: type_,
            rust: crate::search::from_smarts_with_params, python: "from_smarts_with_params", javascript: "fromSmartsWithParams",
            feature: "cap-search", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::SmartsParseParams, default: required }], output: crate::QueryGraph, error: crate::SmartsParseError,
            state: value_returning, operation: none,
            signature: for<'a, 'b> fn(&'a str, &'b crate::SmartsParseParams) -> Result<crate::QueryGraph, crate::SmartsParseError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "search.compile_query", item: callable, owner: module,
            rust: crate::compile_query, python: "compile_query", javascript: "compileQuery",
            feature: "cap-search", status: experimental, kind: module,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }], output: crate::CompiledQuery, error: crate::QueryCompileError,
            state: value_returning, operation: none,
            signature: for<'a> fn(&'a crate::QueryGraph) -> Result<crate::CompiledQuery, crate::QueryCompileError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "search.write_smarts", item: callable, owner: module,
            rust: crate::write_smarts, python: "write_smarts", javascript: "writeSmarts",
            feature: "cap-search", status: experimental, kind: module,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }, { name: params, type: &crate::SmartsWriteParams, default: required }], output: String, error: crate::SmartsWriteError,
            state: value_returning, operation: none,
            signature: for<'a, 'b> fn(&'a crate::QueryGraph, &'b crate::SmartsWriteParams) -> Result<String, crate::SmartsWriteError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "search.write_cx_smarts", item: callable, owner: module,
            rust: crate::write_cx_smarts, python: "write_cx_smarts", javascript: "writeCxSmarts",
            feature: "cap-search", status: experimental, kind: module,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }, { name: params, type: &crate::SmartsWriteParams, default: required }], output: String, error: crate::SmartsWriteError,
            state: value_returning, operation: none,
            signature: for<'a, 'b> fn(&'a crate::QueryGraph, &'b crate::SmartsWriteParams) -> Result<String, crate::SmartsWriteError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "Molecule.substruct_match", item: callable, owner: molecule,
            rust: crate::Molecule::substruct_match, python: "substruct_match", javascript: "substructMatch",
            feature: "cap-search", status: experimental, kind: instance,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }], output: Option<crate::MatchResult>, error: crate::SubstructMatchError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::QueryGraph) -> Result<Option<crate::MatchResult>, crate::SubstructMatchError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "Molecule.substruct_matches", item: callable, owner: molecule,
            rust: crate::Molecule::substruct_matches, python: "substruct_matches", javascript: "substructMatches",
            feature: "cap-search", status: experimental, kind: instance,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }], output: Vec<crate::MatchResult>, error: crate::SubstructMatchError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::QueryGraph) -> Result<Vec<crate::MatchResult>, crate::SubstructMatchError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "Molecule.has_substruct_match", item: callable, owner: molecule,
            rust: crate::Molecule::has_substruct_match, python: "has_substruct_match", javascript: "hasSubstructMatch",
            feature: "cap-search", status: experimental, kind: instance,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }], output: bool, error: crate::SubstructMatchError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::QueryGraph) -> Result<bool, crate::SubstructMatchError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "Molecule.substruct_matches_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::substruct_matches_with_params, python: "substruct_matches_with_params", javascript: "substructMatchesWithParams",
            feature: "cap-search", status: experimental, kind: instance,
            parameters: [{ name: query, type: &crate::QueryGraph, default: required }, { name: params, type: &crate::SubstructMatchParams, default: required }], output: Vec<crate::MatchResult>, error: crate::SubstructMatchError,
            state: read_only, operation: none,
            signature: for<'a, 'b, 'c> fn(&'a crate::Molecule, &'b crate::QueryGraph, &'c crate::SubstructMatchParams) -> Result<Vec<crate::MatchResult>, crate::SubstructMatchError>,
        },
        #[cfg(feature = "cap-search")]
        {
            semantic_id: "Molecule.substruct_matches_compiled", item: callable, owner: molecule,
            rust: crate::Molecule::substruct_matches_compiled, python: "substruct_matches_compiled", javascript: "substructMatchesCompiled",
            feature: "cap-search", status: experimental, kind: instance,
            parameters: [{ name: query, type: &crate::CompiledQuery, default: required }], output: Vec<crate::MatchResult>, error: crate::MatchError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::CompiledQuery) -> Result<Vec<crate::MatchResult>, crate::MatchError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.SmilesWriteParams", item: type, owner: type_,
            rust: crate::SmilesWriteParams, python: "SmilesWriteParams", javascript: "SmilesWriteParams",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.CxSmilesWriteParams", item: type, owner: type_,
            rust: crate::CxSmilesWriteParams, python: "CxSmilesWriteParams", javascript: "CxSmilesWriteParams",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.CxSmilesFields", item: type, owner: type_,
            rust: crate::CxSmilesFields, python: "CxSmilesFields", javascript: "CxSmilesFields",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.CxCoordinateSelection", item: type, owner: type_,
            rust: crate::CxCoordinateSelection, python: "CxCoordinateSelection", javascript: "CxCoordinateSelection",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.RandomSmilesWriteParams", item: type, owner: type_,
            rust: crate::RandomSmilesWriteParams, python: "RandomSmilesWriteParams", javascript: "RandomSmilesWriteParams",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.FragmentSmilesWriteParams", item: type, owner: type_,
            rust: crate::FragmentSmilesWriteParams, python: "FragmentSmilesWriteParams", javascript: "FragmentSmilesWriteParams",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.FragmentCxSmilesWriteParams", item: type, owner: type_,
            rust: crate::FragmentCxSmilesWriteParams, python: "FragmentCxSmilesWriteParams", javascript: "FragmentCxSmilesWriteParams",
            feature: "cap-smiles", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.SmilesWriteError", item: type, owner: type_,
            rust: crate::SmilesWriteError, python: "SmilesWriteError", javascript: "SmilesWriteError",
            feature: "cap-smiles", status: experimental, role: error,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_smiles", item: callable, owner: molecule,
            rust: crate::Molecule::to_smiles, python: "to_smiles", javascript: "toSmiles",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::Molecule) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_smiles_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::to_smiles_with_params, python: "to_smiles_with_params", javascript: "toSmilesWithParams",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::SmilesWriteParams, default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::SmilesWriteParams) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_cx_smiles", item: callable, owner: molecule,
            rust: crate::Molecule::to_cx_smiles, python: "to_cx_smiles", javascript: "toCxSmiles",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::Molecule) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_cx_smiles_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::to_cx_smiles_with_params, python: "to_cx_smiles_with_params", javascript: "toCxSmilesWithParams",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::CxSmilesWriteParams, default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::CxSmilesWriteParams) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_fragment_smiles", item: callable, owner: molecule,
            rust: crate::Molecule::to_fragment_smiles, python: "to_fragment_smiles", javascript: "toFragmentSmiles",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: atoms, type: &[crate::AtomId], default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b [crate::AtomId]) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_fragment_smiles_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::to_fragment_smiles_with_params, python: "to_fragment_smiles_with_params", javascript: "toFragmentSmilesWithParams",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::FragmentSmilesWriteParams, default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::FragmentSmilesWriteParams) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_fragment_cx_smiles", item: callable, owner: molecule,
            rust: crate::Molecule::to_fragment_cx_smiles, python: "to_fragment_cx_smiles", javascript: "toFragmentCxSmiles",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: atoms, type: &[crate::AtomId], default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b [crate::AtomId]) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_fragment_cx_smiles_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::to_fragment_cx_smiles_with_params, python: "to_fragment_cx_smiles_with_params", javascript: "toFragmentCxSmilesWithParams",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::FragmentCxSmilesWriteParams, default: required }],
            output: String, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::FragmentCxSmilesWriteParams) -> Result<String, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_random_smiles", item: callable, owner: molecule,
            rust: crate::Molecule::to_random_smiles, python: "to_random_smiles", javascript: "toRandomSmiles",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: count, type: u32, default: required }, { name: seed, type: u32, default: required }],
            output: Vec<String>, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::Molecule, u32, u32) -> Result<Vec<String>, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.to_random_smiles_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::to_random_smiles_with_params, python: "to_random_smiles_with_params", javascript: "toRandomSmilesWithParams",
            feature: "cap-smiles", status: experimental, kind: instance,
            parameters: [{ name: count, type: u32, default: required }, { name: seed, type: u32, default: required }, { name: params, type: &crate::RandomSmilesWriteParams, default: required }],
            output: Vec<String>, error: crate::SmilesWriteError,
            state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, u32, u32, &'b crate::RandomSmilesWriteParams) -> Result<Vec<String>, crate::SmilesWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSelection", item: type, owner: type_,
            rust: crate::BioSelection, python: "BioSelection", javascript: "BioSelection",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioSelection.from_cid", item: callable, owner: type_,
            rust: crate::BioSelection::from_cid, python: "from_cid", javascript: "fromCid",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: cid, type: &str, default: required }],
            output: crate::BioSelection, error: crate::BioSelectionParseError,
            state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::BioSelection, crate::BioSelectionParseError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.selected_atom_ids", item: callable, owner: type_,
            rust: crate::Protein::selected_atom_ids, python: "selected_atom_ids", javascript: "selectedAtomIds",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: Vec<cosmolkit_bio::BioAtomId>, error: crate::BioSelectionMatchError,
            state: read_only, operation: none,
            signature: fn(&crate::Protein, &crate::BioSelection) -> Result<Vec<cosmolkit_bio::BioAtomId>, crate::BioSelectionMatchError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.selected_atom_ids", item: callable, owner: type_,
            rust: crate::BioStructure::selected_atom_ids, python: "selected_atom_ids", javascript: "selectedAtomIds",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: Vec<cosmolkit_bio::BioAtomId>, error: crate::BioSelectionMatchError,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure, &crate::BioSelection) -> Result<Vec<cosmolkit_bio::BioAtomId>, crate::BioSelectionMatchError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioSelection.to_cid", item: callable, owner: type_,
            rust: crate::BioSelection::to_cid, python: "to_cid", javascript: "toCid",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: String, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioSelection) -> String,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSelectionParseError", item: type, owner: type_,
            rust: crate::BioSelectionParseError, python: "BioSelectionParseError", javascript: "BioSelectionParseError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSelectionMatchError", item: type, owner: type_,
            rust: crate::BioSelectionMatchError, python: "BioSelectionMatchError", javascript: "BioSelectionMatchError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioOperationError", item: type, owner: type_,
            rust: crate::BioOperationError, python: "BioOperationError", javascript: "BioOperationError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.with_translated_coordinates", item: callable, owner: type_,
            rust: crate::BioStructure::with_translated_coordinates, python: "with_translated_coordinates", javascript: "withTranslatedCoordinates",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: offset, type: [f64; 3], default: required }],
            output: crate::BioStructure, error: crate::BioOperationError,
            state: value_returning, operation: "with_translated_coordinates",
            signature: fn(&crate::BioStructure, [f64; 3]) -> Result<crate::BioStructure, crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.translate_", item: callable, owner: type_,
            rust: crate::BioStructure::translate_, python: "translate_", javascript: "translate",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: offset, type: [f64; 3], default: required }],
            output: (), error: crate::BioOperationError,
            state: in_place, operation: "translate_",
            signature: fn(&mut crate::BioStructure, [f64; 3]) -> Result<(), crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.with_translated_coordinates", item: callable, owner: type_,
            rust: crate::Protein::with_translated_coordinates, python: "with_translated_coordinates", javascript: "withTranslatedCoordinates",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: offset, type: [f64; 3], default: required }],
            output: crate::Protein, error: crate::BioOperationError,
            state: value_returning, operation: "with_translated_coordinates",
            signature: fn(&crate::Protein, [f64; 3]) -> Result<crate::Protein, crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.translate_", item: callable, owner: type_,
            rust: crate::Protein::translate_, python: "translate_", javascript: "translate",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: offset, type: [f64; 3], default: required }],
            output: (), error: crate::BioOperationError,
            state: in_place, operation: "translate_",
            signature: fn(&mut crate::Protein, [f64; 3]) -> Result<(), crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.with_selection", item: callable, owner: type_,
            rust: crate::BioStructure::with_selection, python: "with_selection", javascript: "withSelection",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: crate::BioStructure, error: crate::BioOperationError,
            state: value_returning, operation: "with_selection",
            signature: fn(&crate::BioStructure, &crate::BioSelection) -> Result<crate::BioStructure, crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.retain_selection_", item: callable, owner: type_,
            rust: crate::BioStructure::retain_selection_, python: "retain_selection_", javascript: "retainSelection",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: (), error: crate::BioOperationError,
            state: in_place, operation: "retain_selection_",
            signature: fn(&mut crate::BioStructure, &crate::BioSelection) -> Result<(), crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.with_selection", item: callable, owner: type_,
            rust: crate::Protein::with_selection, python: "with_selection", javascript: "withSelection",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: crate::Protein, error: crate::BioOperationError,
            state: value_returning, operation: "with_selection",
            signature: fn(&crate::Protein, &crate::BioSelection) -> Result<crate::Protein, crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.retain_selection_", item: callable, owner: type_,
            rust: crate::Protein::retain_selection_, python: "retain_selection_", javascript: "retainSelection",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: selection, type: &crate::BioSelection, default: required }],
            output: (), error: crate::BioOperationError,
            state: in_place, operation: "retain_selection_",
            signature: fn(&mut crate::Protein, &crate::BioSelection) -> Result<(), crate::BioOperationError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSelectionCopyError", item: type, owner: type_,
            rust: crate::BioSelectionCopyError, python: "BioSelectionCopyError", javascript: "BioSelectionCopyError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSelectionCopyCause", item: type, owner: type_,
            rust: crate::BioSelectionCopyCause, python: "BioSelectionCopyCause", javascript: "BioSelectionCopyCause",
            feature: "cap-bio", status: experimental, role: error,
        },
    #[cfg(feature = "cap-bio")]
    {
        semantic_id: "types.BioRowTraverseError", item: type, owner: type_,
        rust: crate::BioRowTraverseError, python: "BioRowTraverseError", javascript: "BioRowTraverseError",
        feature: "cap-bio", status: experimental, role: error,
    },
    #[cfg(feature = "cap-bio")]
    {
        semantic_id: "types.BioRowModelError", item: type, owner: type_,
        rust: crate::BioRowModelError, python: "BioRowModelError", javascript: "BioRowModelError",
        feature: "cap-bio", status: experimental, role: error,
    },
    #[cfg(feature = "cap-bio")]
    {
        semantic_id: "types.BioRowChainError", item: type, owner: type_,
        rust: crate::BioRowChainError, python: "BioRowChainError", javascript: "BioRowChainError",
        feature: "cap-bio", status: experimental, role: error,

    },
    // Structural readers and explicit Protein projection. Binding names are
        // declarations, not implemented Python/JS adapters; no parity promotion.
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioPdbReadParams", item: type, owner: type_,
            rust: crate::BioPdbReadParams, python: "BioPdbReadParams", javascript: "BioPdbReadParams",
            feature: "cap-bio", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.BioMoleculeParams", item: type, owner: type_,
            rust: crate::BioMoleculeParams, python: "BioMoleculeParams", javascript: "BioMoleculeParams",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.BioMoleculeError", item: type, owner: type_,
            rust: crate::BioMoleculeError, python: "BioMoleculeError", javascript: "BioMoleculeError",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, role: error,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.BioMoleculeConversionError", item: type, owner: type_,
            rust: crate::BioMoleculeConversionError, python: "BioMoleculeConversionError", javascript: "BioMoleculeConversionError",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, role: error,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "BioStructure.to_molecule_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::to_molecule_with_params, python: "to_molecule_with_params", javascript: "toMoleculeWithParams",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, kind: instance,
            parameters: [{name: params, type: &crate::BioMoleculeParams, default: required}],
            output: crate::Molecule, error: crate::BioMoleculeError, state: value_returning, operation: none,
            signature: fn(&crate::BioStructure, &crate::BioMoleculeParams) -> Result<crate::Molecule,crate::BioMoleculeError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "BioStructure.to_molecule", item: callable, owner: type_,
            rust: crate::BioStructure::to_molecule, python: "to_molecule", javascript: "toMolecule",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, kind: instance, parameters: [],
            output: crate::Molecule, error: crate::BioMoleculeError, state: value_returning, operation: none,
            signature: fn(&crate::BioStructure) -> Result<crate::Molecule,crate::BioMoleculeError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "Protein.to_molecule_with_params", item: callable, owner: type_,
            rust: crate::Protein::to_molecule_with_params, python: "to_molecule_with_params", javascript: "toMoleculeWithParams",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, kind: instance,
            parameters: [{name: params, type: &crate::BioMoleculeParams, default: required}],
            output: crate::Molecule, error: crate::BioMoleculeError, state: value_returning, operation: none,
            signature: fn(&crate::Protein, &crate::BioMoleculeParams) -> Result<crate::Molecule,crate::BioMoleculeError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "Protein.to_molecule", item: callable, owner: type_,
            rust: crate::Protein::to_molecule, python: "to_molecule", javascript: "toMolecule",
            feature: "cap-io", requires: ["cap-bio"], status: experimental, kind: instance, parameters: [],
            output: crate::Molecule, error: crate::BioMoleculeError, state: value_returning, operation: none,
            signature: fn(&crate::Protein) -> Result<crate::Molecule,crate::BioMoleculeError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioPdbReadError", item: type, owner: type_,
            rust: crate::BioPdbReadError, python: "BioPdbReadError", javascript: "BioPdbReadError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioPdbReadStage", item: type, owner: type_,
            rust: crate::BioPdbReadStage, python: "BioPdbReadStage", javascript: "BioPdbReadStage",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioMmcifReadError", item: type, owner: type_,
            rust: crate::BioMmcifReadError, python: "BioMmcifReadError", javascript: "BioMmcifReadError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioMmcifReadStage", item: type, owner: type_,
            rust: crate::BioMmcifReadStage, python: "BioMmcifReadStage", javascript: "BioMmcifReadStage",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioReadParams", item: type, owner: type_,
            rust: crate::BioReadParams, python: "BioReadParams", javascript: "BioReadParams",
            feature: "cap-bio", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioReadError", item: type, owner: type_,
            rust: crate::BioReadError, python: "BioReadError", javascript: "BioReadError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinReadError", item: type, owner: type_,
            rust: crate::ProteinReadError, python: "ProteinReadError", javascript: "ProteinReadError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioCrystalInfo.space_group_number", item: callable, owner: type_,
            rust: crate::BioCrystalInfo::space_group_number, python: "space_group_number", javascript: "spaceGroupNumber",
            feature: "cap-bio", status: experimental, kind: instance, parameters: [],
            output: Option<i32>, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioCrystalInfo) -> Option<i32>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioTransform.approx", item: callable, owner: type_,
            rust: crate::BioTransform::approx, python: "approx", javascript: "approx",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: other, type: &crate::BioTransform, default: required }, { name: epsilon, type: f64, default: required }],
            output: bool, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioTransform, &crate::BioTransform, f64) -> bool,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioMmcifWriteParams", item: type, owner: type_,
            rust: crate::BioMmcifWriteParams, python: "BioMmcifWriteParams", javascript: "BioMmcifWriteParams",
            feature: "cap-bio", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioMmcifWriteError", item: type, owner: type_,
            rust: crate::BioMmcifWriteError, python: "BioMmcifWriteError", javascript: "BioMmcifWriteError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.to_mmcif_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::to_mmcif_with_params, python: "to_mmcif_with_params", javascript: "toMmcifWithParams",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::BioMmcifWriteParams, default: required }],
            output: String, error: crate::BioMmcifWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &crate::BioMmcifWriteParams) -> Result<String, crate::BioMmcifWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.to_mmcif", item: callable, owner: type_,
            rust: crate::BioStructure::to_mmcif, python: "to_mmcif", javascript: "toMmcif",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: String, error: crate::BioMmcifWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> Result<String, crate::BioMmcifWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.write_mmcif_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::write_mmcif_with_params, python: "write_mmcif_with_params", javascript: "writeMmcifWithParams",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: path, type: &std::path::Path, default: required }, { name: params, type: &crate::BioMmcifWriteParams, default: required }],
            output: (), error: crate::BioMmcifWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &std::path::Path, &crate::BioMmcifWriteParams) -> Result<(), crate::BioMmcifWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.write_mmcif", item: callable, owner: type_,
            rust: crate::BioStructure::write_mmcif, python: "write_mmcif", javascript: "writeMmcif",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: path, type: &std::path::Path, default: required }],
            output: (), error: crate::BioMmcifWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &std::path::Path) -> Result<(), crate::BioMmcifWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioPdbWriteParams", item: type, owner: type_,
            rust: crate::BioPdbWriteParams, python: "BioPdbWriteParams", javascript: "BioPdbWriteParams",
            feature: "cap-bio", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioPdbWriteError", item: type, owner: type_,
            rust: crate::BioPdbWriteError, python: "BioPdbWriteError", javascript: "BioPdbWriteError",
            feature: "cap-bio", status: experimental, role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.to_pdb_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::to_pdb_with_params, python: "to_pdb_with_params", javascript: "toPdbWithParams",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::BioPdbWriteParams, default: required }],
            output: String, error: crate::BioPdbWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &crate::BioPdbWriteParams) -> Result<String, crate::BioPdbWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.to_pdb", item: callable, owner: type_,
            rust: crate::BioStructure::to_pdb, python: "to_pdb", javascript: "toPdb",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: String, error: crate::BioPdbWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> Result<String, crate::BioPdbWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.write_pdb_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::write_pdb_with_params, python: "write_pdb_with_params", javascript: "writePdbWithParams",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: path, type: &std::path::Path, default: required }, { name: params, type: &crate::BioPdbWriteParams, default: required }],
            output: (), error: crate::BioPdbWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &std::path::Path, &crate::BioPdbWriteParams) -> Result<(), crate::BioPdbWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.write_pdb", item: callable, owner: type_,
            rust: crate::BioStructure::write_pdb, python: "write_pdb", javascript: "writePdb",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: path, type: &std::path::Path, default: required }],
            output: (), error: crate::BioPdbWriteError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, &std::path::Path) -> Result<(), crate::BioPdbWriteError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_text_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::from_text_with_params, python: "from_text_with_params", javascript: "fromTextWithParams",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::BioReadParams, default: required }],
            output: crate::BioStructure, error: crate::BioReadError, state: value_returning, operation: none,
            signature: fn(&str, &crate::BioReadParams) -> Result<crate::BioStructure, crate::BioReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_text", item: callable, owner: type_,
            rust: crate::BioStructure::from_text, python: "from_text", javascript: "fromText",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::BioStructure, error: crate::BioReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::BioStructure, crate::BioReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.read_with_format", item: callable, owner: type_,
            rust: crate::BioStructure::read_with_format, python: "read_with_format", javascript: "readWithFormat",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: path, type: &std::path::Path, default: required }, { name: format, type: crate::BioCoordinateFormat, default: required }],
            output: crate::BioStructure, error: crate::BioReadError, state: value_returning, operation: none,
            signature: fn(&std::path::Path, crate::BioCoordinateFormat) -> Result<crate::BioStructure, crate::BioReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.read", item: callable, owner: type_,
            rust: crate::BioStructure::read, python: "read", javascript: "read",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: path, type: &std::path::Path, default: required }],
            output: crate::BioStructure, error: crate::BioReadError, state: value_returning, operation: none,
            signature: fn(&std::path::Path) -> Result<crate::BioStructure, crate::BioReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.from_text_with_params", item: callable, owner: type_,
            rust: crate::Protein::from_text_with_params, python: "from_text_with_params", javascript: "fromTextWithParams",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::BioReadParams, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&str, &crate::BioReadParams) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.from_text", item: callable, owner: type_,
            rust: crate::Protein::from_text, python: "from_text", javascript: "fromText",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.read_with_format", item: callable, owner: type_,
            rust: crate::Protein::read_with_format, python: "read_with_format", javascript: "readWithFormat",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: path, type: &std::path::Path, default: required }, { name: format, type: crate::BioCoordinateFormat, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&std::path::Path, crate::BioCoordinateFormat) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.read", item: callable, owner: type_,
            rust: crate::Protein::read, python: "read", javascript: "read",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: path, type: &std::path::Path, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&std::path::Path) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_pdb", item: callable, owner: type_,
            rust: crate::BioStructure::from_pdb, python: "from_pdb", javascript: "fromPdb",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::BioStructure, error: crate::BioPdbReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::BioStructure, crate::BioPdbReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_pdb_with_params", item: callable, owner: type_,
            rust: crate::BioStructure::from_pdb_with_params, python: "from_pdb_with_params", javascript: "fromPdbWithParams",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::BioPdbReadParams, default: required }],
            output: crate::BioStructure, error: crate::BioPdbReadError, state: value_returning, operation: none,
            signature: fn(&str, &crate::BioPdbReadParams) -> Result<crate::BioStructure, crate::BioPdbReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_mmcif", item: callable, owner: type_,
            rust: crate::BioStructure::from_mmcif, python: "from_mmcif", javascript: "fromMmcif",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::BioStructure, error: crate::BioMmcifReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::BioStructure, crate::BioMmcifReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.from_pdb", item: callable, owner: type_,
            rust: crate::Protein::from_pdb, python: "from_pdb", javascript: "fromPdb",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.from_pdb_with_params", item: callable, owner: type_,
            rust: crate::Protein::from_pdb_with_params, python: "from_pdb_with_params", javascript: "fromPdbWithParams",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }, { name: params, type: &crate::BioPdbReadParams, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&str, &crate::BioPdbReadParams) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.from_mmcif", item: callable, owner: type_,
            rust: crate::Protein::from_mmcif, python: "from_mmcif", javascript: "fromMmcif",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: text, type: &str, default: required }],
            output: crate::Protein, error: crate::ProteinReadError, state: value_returning, operation: none,
            signature: fn(&str) -> Result<crate::Protein, crate::ProteinReadError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioPdbReadError.stage", item: callable, owner: type_,
            rust: crate::BioPdbReadError::stage, python: "stage", javascript: "stage",
            feature: "cap-bio", status: experimental, kind: instance, parameters: [],
            output: crate::BioPdbReadStage, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioPdbReadError) -> crate::BioPdbReadStage,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioPdbReadError.line_number", item: callable, owner: type_,
            rust: crate::BioPdbReadError::line_number, python: "line_number", javascript: "lineNumber",
            feature: "cap-bio", status: experimental, kind: instance, parameters: [],
            output: Option<i32>, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioPdbReadError) -> Option<i32>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioPdbReadError.record_tag", item: callable, owner: type_,
            rust: crate::BioPdbReadError::record_tag, python: "record_tag", javascript: "recordTag",
            feature: "cap-bio", status: experimental, kind: instance, parameters: [],
            output: Option<[u8; 4]>, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioPdbReadError) -> Option<[u8; 4]>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioMmcifReadError.stage", item: callable, owner: type_,
            rust: crate::BioMmcifReadError::stage, python: "stage", javascript: "stage",
            feature: "cap-bio", status: experimental, kind: instance, parameters: [],
            output: crate::BioMmcifReadStage, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioMmcifReadError) -> crate::BioMmcifReadStage,
        },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomPairParams", item: type, owner: type_, rust: crate::AtomPairParams, python: "AtomPairParams", javascript: "AtomPairParams", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomPairFingerprintParams", item: type, owner: type_, rust: crate::AtomPairFingerprintParams, python: "AtomPairFingerprintParams", javascript: "AtomPairFingerprintParams", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomPairAtomInvariantsGenerator", item: type, owner: type_, rust: crate::AtomPairAtomInvariantsGenerator, python: "AtomPairAtomInvariantsGenerator", javascript: "AtomPairAtomInvariantsGenerator", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomPairReadError", item: type, owner: type_, rust: crate::AtomPairReadError, python: "AtomPairReadError", javascript: "AtomPairReadError", feature: "cap-fingerprints", status: experimental, role: error, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_fingerprint, python: "atom_pair_fingerprint", javascript: "atomPairFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::Fingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::Fingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_fingerprint_with_params, python: "atom_pair_fingerprint_with_params", javascript: "atomPairFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::AtomPairFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::Fingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::AtomPairFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::Fingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_sparse_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_sparse_fingerprint, python: "atom_pair_sparse_fingerprint", javascript: "atomPairSparseFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseBitFingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseBitFingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_sparse_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_sparse_fingerprint_with_params, python: "atom_pair_sparse_fingerprint_with_params", javascript: "atomPairSparseFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::AtomPairFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseBitFingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::AtomPairFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseBitFingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_count_fingerprint, python: "atom_pair_count_fingerprint", javascript: "atomPairCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint32, error: crate::AtomPairReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseCountFingerprint32,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_count_fingerprint_with_params, python: "atom_pair_count_fingerprint_with_params", javascript: "atomPairCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::AtomPairFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseCountFingerprint32, error: crate::AtomPairReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::AtomPairFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseCountFingerprint32,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_sparse_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_sparse_count_fingerprint, python: "atom_pair_sparse_count_fingerprint", javascript: "atomPairSparseCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseCountFingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.atom_pair_sparse_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::atom_pair_sparse_count_fingerprint_with_params, python: "atom_pair_sparse_count_fingerprint_with_params", javascript: "atomPairSparseCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::AtomPairFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseCountFingerprint, error: crate::AtomPairReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::AtomPairFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseCountFingerprint,crate::AtomPairReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairAtomInvariantsGenerator.info_string", item: callable, owner: type_, rust: crate::AtomPairAtomInvariantsGenerator::info_string, python: "info_string", javascript: "infoString", feature: "cap-fingerprints", status: experimental, kind: instance, parameters: [], output: String, error: none, state: read_only, operation: none, signature: fn(&crate::AtomPairAtomInvariantsGenerator) -> String, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomPairAtomInvariantsGenerator.to_json", item: callable, owner: type_, rust: crate::AtomPairAtomInvariantsGenerator::to_json, python: "to_json", javascript: "toJson", feature: "cap-fingerprints", status: experimental, kind: instance, parameters: [], output: String, error: none, state: read_only, operation: none, signature: fn(&crate::AtomPairAtomInvariantsGenerator) -> String, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.AtomCodeExplanation", item: type, owner: type_, rust: crate::AtomCodeExplanation, python: "AtomCodeExplanation", javascript: "AtomCodeExplanation", feature: "cap-fingerprints", status: experimental, role: result, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "errors.AtomCodeExplanationError", item: type, owner: type_, rust: crate::AtomCodeExplanationError, python: "AtomCodeExplanationError", javascript: "AtomCodeExplanationError", feature: "cap-fingerprints", status: experimental, role: error, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.TopologicalTorsionParams", item: type, owner: type_, rust: crate::TopologicalTorsionParams, python: "TopologicalTorsionParams", javascript: "TopologicalTorsionParams", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.TopologicalTorsionFingerprintParams", item: type, owner: type_, rust: crate::TopologicalTorsionFingerprintParams, python: "TopologicalTorsionFingerprintParams", javascript: "TopologicalTorsionFingerprintParams", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.TopologicalTorsionReadError", item: type, owner: type_, rust: crate::TopologicalTorsionReadError, python: "TopologicalTorsionReadError", javascript: "TopologicalTorsionReadError", feature: "cap-fingerprints", status: experimental, role: error, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_fingerprint, python: "topological_torsion_fingerprint", javascript: "topologicalTorsionFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::Fingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::Fingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_fingerprint_with_params, python: "topological_torsion_fingerprint_with_params", javascript: "topologicalTorsionFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::TopologicalTorsionFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::Fingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::TopologicalTorsionFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::Fingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_sparse_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_sparse_fingerprint, python: "topological_torsion_sparse_fingerprint", javascript: "topologicalTorsionSparseFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseBitFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseBitFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_sparse_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_sparse_fingerprint_with_params, python: "topological_torsion_sparse_fingerprint_with_params", javascript: "topologicalTorsionSparseFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::TopologicalTorsionFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseBitFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::TopologicalTorsionFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseBitFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_count_fingerprint, python: "topological_torsion_count_fingerprint", javascript: "topologicalTorsionCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint32, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseCountFingerprint32,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_count_fingerprint_with_params, python: "topological_torsion_count_fingerprint_with_params", javascript: "topologicalTorsionCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::TopologicalTorsionFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseCountFingerprint32, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::TopologicalTorsionFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseCountFingerprint32,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_sparse_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_sparse_count_fingerprint, python: "topological_torsion_sparse_count_fingerprint", javascript: "topologicalTorsionSparseCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_sparse_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_sparse_count_fingerprint_with_params, python: "topological_torsion_sparse_count_fingerprint_with_params", javascript: "topologicalTorsionSparseCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::TopologicalTorsionFingerprintParams, default: required }, { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required }], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b,'c> fn(&'a crate::Molecule, &'b crate::TopologicalTorsionFingerprintParams, Option<&'c mut crate::FingerprintAdditionalOutput>) -> Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.FingerprintJsonError",item:type,owner:type_,rust:crate::FingerprintJsonError,python:"FingerprintJsonError",javascript:"FingerprintJsonError",feature:"cap-fingerprints",status:experimental,role:error},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"MorganParams.info_string",item:callable,owner:type_,rust:crate::MorganParams::info_string,python:"info_string",javascript:"infoString",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MorganParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"MorganParams.to_json",item:callable,owner:type_,rust:crate::MorganParams::to_json,python:"to_json",javascript:"toJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MorganParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"MorganParams.with_json",item:callable,owner:type_,rust:crate::MorganParams::with_json,python:"with_json",javascript:"withJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:json,type:&str,default:required}],output:crate::MorganParams,error:crate::FingerprintJsonError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::MorganParams,&'b str)->Result<crate::MorganParams,crate::FingerprintJsonError>},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"AtomPairParams.info_string",item:callable,owner:type_,rust:crate::AtomPairParams::info_string,python:"info_string",javascript:"infoString",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::AtomPairParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"AtomPairParams.to_json",item:callable,owner:type_,rust:crate::AtomPairParams::to_json,python:"to_json",javascript:"toJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::AtomPairParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"AtomPairParams.with_json",item:callable,owner:type_,rust:crate::AtomPairParams::with_json,python:"with_json",javascript:"withJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:json,type:&str,default:required}],output:crate::AtomPairParams,error:crate::FingerprintJsonError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::AtomPairParams,&'b str)->Result<crate::AtomPairParams,crate::FingerprintJsonError>},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"TopologicalTorsionParams.info_string",item:callable,owner:type_,rust:crate::TopologicalTorsionParams::info_string,python:"info_string",javascript:"infoString",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::TopologicalTorsionParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"TopologicalTorsionParams.to_json",item:callable,owner:type_,rust:crate::TopologicalTorsionParams::to_json,python:"to_json",javascript:"toJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:String,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::TopologicalTorsionParams)->String},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"TopologicalTorsionParams.with_json",item:callable,owner:type_,rust:crate::TopologicalTorsionParams::with_json,python:"with_json",javascript:"withJson",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:json,type:&str,default:required}],output:crate::TopologicalTorsionParams,error:crate::FingerprintJsonError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::TopologicalTorsionParams,&'b str)->Result<crate::TopologicalTorsionParams,crate::FingerprintJsonError>},
        // FP-values public boundary. Source: pinned RDKit SparseIntVect.h
        // and Wrap/SparseIntVect.cpp:135-185. Canonical names replace get_*
        // and *_i implementation spellings. Python/JS projections are declared,
        // not implemented by this Rust-only registration. No helper exports.
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.MorganParams", item: type, owner: type_,
            rust: crate::MorganParams, python: "MorganParams", javascript: "MorganParams",
            feature: "cap-fingerprints", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.MorganInvariants", item: type, owner: type_,
            rust: crate::MorganInvariants, python: "MorganInvariants", javascript: "MorganInvariants",
            feature: "cap-fingerprints", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.MorganFingerprintParams", item: type, owner: type_,
            rust: crate::MorganFingerprintParams, python: "MorganFingerprintParams", javascript: "MorganFingerprintParams",
            feature: "cap-fingerprints", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.MorganReadError", item: type, owner: type_,
            rust: crate::MorganReadError, python: "MorganReadError", javascript: "MorganReadError",
            feature: "cap-fingerprints", status: experimental, role: error,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.FingerprintAdditionalOutput", item: type, owner: type_,
            rust: crate::FingerprintAdditionalOutput, python: "FingerprintAdditionalOutput", javascript: "FingerprintAdditionalOutput",
            feature: "cap-fingerprints", status: experimental, role: value,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.Fingerprint", item: type, owner: type_,
            rust: crate::Fingerprint, python: "Fingerprint", javascript: "Fingerprint",
            feature: "cap-fingerprints", status: experimental, role: value,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.SparseBitFingerprint", item: type, owner: type_,
            rust: crate::SparseBitFingerprint, python: "SparseBitFingerprint", javascript: "SparseBitFingerprint",
            feature: "cap-fingerprints", status: experimental, role: value,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.SparseCountFingerprint", item: type, owner: type_,
            rust: crate::SparseCountFingerprint, python: "SparseCountFingerprint", javascript: "SparseCountFingerprint",
            feature: "cap-fingerprints", status: experimental, role: value,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.SparseCountFingerprint32", item: type, owner: type_,
            rust: crate::SparseCountFingerprint32, python: "SparseCountFingerprint32", javascript: "SparseCountFingerprint32",
            feature: "cap-fingerprints", status: experimental, role: value,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "types.FingerprintError", item: type, owner: type_,
            rust: crate::FingerprintError, python: "FingerprintError", javascript: "FingerprintError",
            feature: "cap-fingerprints", status: experimental, role: error,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_sparse_count_fingerprint",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_sparse_count_fingerprint,
            python: "morgan_sparse_count_fingerprint",
            javascript: "morganSparseCountFingerprint",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [],
            output: crate::SparseCountFingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::Molecule,
            ) -> Result<crate::SparseCountFingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_sparse_count_fingerprint_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_sparse_count_fingerprint_with_params,
            python: "morgan_sparse_count_fingerprint_with_params",
            javascript: "morganSparseCountFingerprintWithParams",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [
                { name: params, type: &crate::MorganFingerprintParams, default: required },
                { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required },
            ],
            output: crate::SparseCountFingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b, 'c> fn(
                &'a crate::Molecule,
                &'b crate::MorganFingerprintParams,
                Option<&'c mut crate::FingerprintAdditionalOutput>,
            ) -> Result<crate::SparseCountFingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_sparse_fingerprint",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_sparse_fingerprint,
            python: "morgan_sparse_fingerprint",
            javascript: "morganSparseFingerprint",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [],
            output: crate::SparseBitFingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::Molecule,
            ) -> Result<crate::SparseBitFingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_sparse_fingerprint_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_sparse_fingerprint_with_params,
            python: "morgan_sparse_fingerprint_with_params",
            javascript: "morganSparseFingerprintWithParams",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [
                { name: params, type: &crate::MorganFingerprintParams, default: required },
                { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required },
            ],
            output: crate::SparseBitFingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b, 'c> fn(
                &'a crate::Molecule,
                &'b crate::MorganFingerprintParams,
                Option<&'c mut crate::FingerprintAdditionalOutput>,
            ) -> Result<crate::SparseBitFingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_count_fingerprint",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_count_fingerprint,
            python: "morgan_count_fingerprint",
            javascript: "morganCountFingerprint",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [],
            output: crate::SparseCountFingerprint32,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::Molecule,
            ) -> Result<crate::SparseCountFingerprint32, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_count_fingerprint_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_count_fingerprint_with_params,
            python: "morgan_count_fingerprint_with_params",
            javascript: "morganCountFingerprintWithParams",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [
                { name: params, type: &crate::MorganFingerprintParams, default: required },
                { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required },
            ],
            output: crate::SparseCountFingerprint32,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b, 'c> fn(
                &'a crate::Molecule,
                &'b crate::MorganFingerprintParams,
                Option<&'c mut crate::FingerprintAdditionalOutput>,
            ) -> Result<crate::SparseCountFingerprint32, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_fingerprint",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_fingerprint,
            python: "morgan_fingerprint",
            javascript: "morganFingerprint",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [],
            output: crate::Fingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::Molecule,
            ) -> Result<crate::Fingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.morgan_fingerprint_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::morgan_fingerprint_with_params,
            python: "morgan_fingerprint_with_params",
            javascript: "morganFingerprintWithParams",
            feature: "cap-fingerprints", status: experimental,
            kind: instance, receiver: shared,
            parameters: [
                { name: params, type: &crate::MorganFingerprintParams, default: required },
                { name: additional_output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: required },
            ],
            output: crate::Fingerprint,
            error: crate::MorganReadError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b, 'c> fn(
                &'a crate::Molecule,
                &'b crate::MorganFingerprintParams,
                Option<&'c mut crate::FingerprintAdditionalOutput>,
            ) -> Result<crate::Fingerprint, crate::MorganReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.new", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::new, python: "new", javascript: "new",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [], output: crate::FingerprintAdditionalOutput, error: none,
            state: value_returning, operation: none,
            signature: fn() -> crate::FingerprintAdditionalOutput,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.default", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::default, python: "default", javascript: "default",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [],
            output: crate::FingerprintAdditionalOutput, error: none,
            state: value_returning, operation: none,
            signature: fn() -> crate::FingerprintAdditionalOutput,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.allocate_atom_counts", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::allocate_atom_counts, python: "allocate_atom_counts", javascript: "allocateAtomCounts",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [], output: (), error: none,
            state: in_place, operation: none,
            signature: fn(&mut crate::FingerprintAdditionalOutput),
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.allocate_atom_to_bits", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::allocate_atom_to_bits, python: "allocate_atom_to_bits", javascript: "allocateAtomToBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [], output: (), error: none,
            state: in_place, operation: none,
            signature: fn(&mut crate::FingerprintAdditionalOutput),
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.allocate_bit_info_map", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::allocate_bit_info_map, python: "allocate_bit_info_map", javascript: "allocateBitInfoMap",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [], output: (), error: none,
            state: in_place, operation: none,
            signature: fn(&mut crate::FingerprintAdditionalOutput),
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.allocate_bit_paths", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::allocate_bit_paths, python: "allocate_bit_paths", javascript: "allocateBitPaths",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [], output: (), error: none,
            state: in_place, operation: none,
            signature: fn(&mut crate::FingerprintAdditionalOutput),
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.allocate_atoms_per_bit", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::allocate_atoms_per_bit, python: "allocate_atoms_per_bit", javascript: "allocateAtomsPerBit",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [], output: (), error: none,
            state: in_place, operation: none,
            signature: fn(&mut crate::FingerprintAdditionalOutput),
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.atom_counts", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::atom_counts, python: "atom_counts", javascript: "atomCounts",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Option<&[u32]>, error: none,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::FingerprintAdditionalOutput) -> Option<&'a [u32]>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.atom_to_bits", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::atom_to_bits, python: "atom_to_bits", javascript: "atomToBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Option<&[Vec<u64>]>, error: none,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::FingerprintAdditionalOutput) -> Option<&'a [Vec<u64>]>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.bit_info_map", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::bit_info_map, python: "bit_info_map", javascript: "bitInfoMap",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Option<&std::collections::BTreeMap<u64, Vec<(u32, u32)>>>, error: none,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::FingerprintAdditionalOutput) -> Option<&'a std::collections::BTreeMap<u64, Vec<(u32, u32)>>>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.bit_paths", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::bit_paths, python: "bit_paths", javascript: "bitPaths",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Option<&std::collections::BTreeMap<u64, Vec<Vec<i32>>>>, error: none,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::FingerprintAdditionalOutput) -> Option<&'a std::collections::BTreeMap<u64, Vec<Vec<i32>>>>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "FingerprintAdditionalOutput.atoms_per_bit", item: callable, owner: type_,
            rust: crate::FingerprintAdditionalOutput::atoms_per_bit, python: "atoms_per_bit", javascript: "atomsPerBit",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Option<&std::collections::BTreeMap<u64, Vec<Vec<i32>>>>, error: none,
            state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::FingerprintAdditionalOutput) -> Option<&'a std::collections::BTreeMap<u64, Vec<Vec<i32>>>>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Fingerprint.n_bits", item: callable, owner: type_,
            rust: crate::Fingerprint::n_bits, python: "n_bits", javascript: "nBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: u32, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::Fingerprint) -> u32,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Fingerprint.on_bits", item: callable, owner: type_,
            rust: crate::Fingerprint::on_bits, python: "on_bits", javascript: "onBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Vec<u32>, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::Fingerprint) -> Vec<u32>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseBitFingerprint.n_bits", item: callable, owner: type_,
            rust: crate::SparseBitFingerprint::n_bits, python: "n_bits", javascript: "nBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: u32, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseBitFingerprint) -> u32,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseBitFingerprint.on_bits", item: callable, owner: type_,
            rust: crate::SparseBitFingerprint::on_bits, python: "on_bits", javascript: "onBits",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: Vec<i32>, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseBitFingerprint) -> Vec<i32>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.new", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::new, python: "new", javascript: "new",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [{ name: length, type: u64, default: required }],
            output: crate::SparseCountFingerprint, error: none,
            state: value_returning, operation: none,
            signature: fn(u64) -> crate::SparseCountFingerprint,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.length", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::length, python: "length", javascript: "length",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [],
            output: u64, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint) -> u64,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::value, python: "value", javascript: "value",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: index, type: u64, default: required }],
            output: i32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, u64) -> Result<i32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.set_value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::set_value, python: "set_value", javascript: "setValue",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: index, type: u64, default: required }, { name: value, type: i32, default: required }],
            output: (), error: crate::FingerprintError,
            state: in_place, operation: none,
            signature: fn(&mut crate::SparseCountFingerprint, u64, i32) -> Result<(), crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.nonzero_elements", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::nonzero_elements, python: "nonzero_elements", javascript: "nonzeroElements",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [],
            output: &std::collections::BTreeMap<u64, i32>, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint) -> &std::collections::BTreeMap<u64, i32>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.total_value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::total_value, python: "total_value", javascript: "totalValue",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: use_abs, type: bool, default: false }],
            output: i32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, bool) -> Result<i32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.fuzzy_and", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::fuzzy_and, python: "fuzzy_and", javascript: "fuzzyAnd",
            feature: "cap-fingerprints", status: parity("RDKit"), kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, &crate::SparseCountFingerprint) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.fuzzy_or", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::fuzzy_or, python: "fuzzy_or", javascript: "fuzzyOr",
            feature: "cap-fingerprints", status: parity("RDKit"), kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, &crate::SparseCountFingerprint) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_added", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_added, python: "with_added", javascript: "withAdded",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, &crate::SparseCountFingerprint) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_subtracted", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_subtracted, python: "with_subtracted", javascript: "withSubtracted",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, &crate::SparseCountFingerprint) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_added_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_added_scalar, python: "with_added_scalar", javascript: "withAddedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, i32) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_subtracted_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_subtracted_scalar, python: "with_subtracted_scalar", javascript: "withSubtractedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, i32) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_multiplied_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_multiplied_scalar, python: "with_multiplied_scalar", javascript: "withMultipliedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, i32) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint.with_divided_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint::with_divided_scalar, python: "with_divided_scalar", javascript: "withDividedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint, i32) -> Result<crate::SparseCountFingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.new", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::new, python: "new", javascript: "new",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [{ name: length, type: u32, default: required }],
            output: crate::SparseCountFingerprint32, error: none,
            state: value_returning, operation: none,
            signature: fn(u32) -> crate::SparseCountFingerprint32,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.length", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::length, python: "length", javascript: "length",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [],
            output: u32, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32) -> u32,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::value, python: "value", javascript: "value",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: index, type: u32, default: required }],
            output: i32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, u32) -> Result<i32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.set_value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::set_value, python: "set_value", javascript: "setValue",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: index, type: u32, default: required }, { name: value, type: i32, default: required }],
            output: (), error: crate::FingerprintError,
            state: in_place, operation: none,
            signature: fn(&mut crate::SparseCountFingerprint32, u32, i32) -> Result<(), crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.nonzero_elements", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::nonzero_elements, python: "nonzero_elements", javascript: "nonzeroElements",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [],
            output: &std::collections::BTreeMap<u32, i32>, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32) -> &std::collections::BTreeMap<u32, i32>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.total_value", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::total_value, python: "total_value", javascript: "totalValue",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: use_abs, type: bool, default: false }],
            output: i32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, bool) -> Result<i32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.fuzzy_and", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::fuzzy_and, python: "fuzzy_and", javascript: "fuzzyAnd",
            feature: "cap-fingerprints", status: parity("RDKit"), kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, &crate::SparseCountFingerprint32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.fuzzy_or", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::fuzzy_or, python: "fuzzy_or", javascript: "fuzzyOr",
            feature: "cap-fingerprints", status: parity("RDKit"), kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, &crate::SparseCountFingerprint32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_added", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_added, python: "with_added", javascript: "withAdded",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, &crate::SparseCountFingerprint32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_subtracted", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_subtracted, python: "with_subtracted", javascript: "withSubtracted",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: other, type: &crate::SparseCountFingerprint32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, &crate::SparseCountFingerprint32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_added_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_added_scalar, python: "with_added_scalar", javascript: "withAddedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, i32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_subtracted_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_subtracted_scalar, python: "with_subtracted_scalar", javascript: "withSubtractedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, i32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_multiplied_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_multiplied_scalar, python: "with_multiplied_scalar", javascript: "withMultipliedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, i32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "SparseCountFingerprint32.with_divided_scalar", item: callable, owner: type_,
            rust: crate::SparseCountFingerprint32::with_divided_scalar, python: "with_divided_scalar", javascript: "withDividedScalar",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: value, type: i32, default: required }],
            output: crate::SparseCountFingerprint32, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::SparseCountFingerprint32, i32) -> Result<crate::SparseCountFingerprint32, crate::FingerprintError>,
        },

                {
            semantic_id:"types.AtomMetadata",item:type,owner:type_,
            rust:crate::AtomMetadata,python:"AtomMetadata",javascript:"AtomMetadata",
            feature:"metadata",status:experimental,role:result,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id:"Molecule.atom_metadata",item:callable,owner:molecule,
            rust:crate::Molecule::atom_metadata,python:"atom_metadata",javascript:"atomMetadata",
            feature:"cap-valence",status:experimental,kind:instance,
            parameters:[{name:recalculate,type:bool,default:boolean(true)}],
            output:Vec<crate::AtomMetadata>,error:crate::ValenceError,state:read_only,operation:none,
            signature:fn(&crate::Molecule,bool)->Result<Vec<crate::AtomMetadata>,crate::ValenceError>,
        },
        // RUN-state owns exactly this public value entry. Its Arc-backed
        // MoleculeState and derived-cache authority stay private; construction
        // and read callables are registered by RUN-builder and RUN-read.
        {
            semantic_id: "types.Molecule",
            item: type,
            owner: type_,
            rust: crate::Molecule,
            python: "Molecule",
            javascript: "Molecule",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.MoleculeBuilder",
            item: type,
            owner: type_,
            rust: crate::MoleculeBuilder,
            python: "MoleculeBuilder",
            javascript: "MoleculeBuilder",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "Molecule.new",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [],
            output: crate::Molecule,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn() -> crate::Molecule,
        },
        {
            semantic_id: "Molecule.from_parts",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::from_parts,
            python: "from_parts",
            javascript: "fromParts",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [
                { name: topology, type: crate::TopologyBlock, default: required },
                { name: coordinates, type: crate::CoordinateBlock, default: required },
                { name: properties, type: crate::MoleculeProperties, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: none,
            signature: fn(crate::TopologyBlock, crate::CoordinateBlock, crate::MoleculeProperties)
                -> Result<crate::Molecule, crate::OperationError>,
        },
        {
            semantic_id: "Molecule.to_builder",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::to_builder,
            python: "to_builder",
            javascript: "toBuilder",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::MoleculeBuilder,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "Molecule.num_atoms",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_atoms,
            python: "num_atoms",
            javascript: "numAtoms",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> usize,
        },
        {
            semantic_id: "Molecule.num_bonds",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_bonds,
            python: "num_bonds",
            javascript: "numBonds",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> usize,
        },
        {
            semantic_id: "Molecule.atoms",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::atoms,
            python: "atoms",
            javascript: "atoms",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::Atom],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> &[crate::Atom],
        },
        {
            semantic_id: "Molecule.bonds",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::bonds,
            python: "bonds",
            javascript: "bonds",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::Bond],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> &[crate::Bond],
        },
        {
            semantic_id: "Molecule.atom",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::atom,
            python: "atom",
            javascript: "atom",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: atom_id, type: crate::AtomId, default: required }],
            output: Option<&crate::Atom>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule, crate::AtomId) -> Option<&crate::Atom>,
        },
        {
            semantic_id: "Molecule.bond",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::bond,
            python: "bond",
            javascript: "bond",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: bond_id, type: crate::BondId, default: required }],
            output: Option<&crate::Bond>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule, crate::BondId) -> Option<&crate::Bond>,
        },
        {
            semantic_id: "Molecule.topology",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::topology,
            python: "topology",
            javascript: "topology",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::TopologyBlock,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> &crate::TopologyBlock,
        },
        {
            semantic_id: "Molecule.coordinates_2d",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::coordinates_2d,
            python: "coordinates_2d",
            javascript: "coordinates2d",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: Option<&[[f64; 2]]>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Option<&[[f64; 2]]>,
        },
        {
            semantic_id: "Molecule.conformers_3d",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::conformers_3d,
            python: "conformers_3d",
            javascript: "conformers3d",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::Conformer3D],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> &[crate::Conformer3D],
        },
        {
            semantic_id: "Molecule.properties",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::properties,
            python: "properties",
            javascript: "properties",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::MoleculeProperties,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> &crate::MoleculeProperties,
        },
        {
            semantic_id: "Molecule.property",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::property,
            python: "property",
            javascript: "property",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: key, type: &str, default: required }],
            output: Option<&str>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b str) -> Option<&'a str>,
        },
        {
            semantic_id: "MoleculeBuilder.new",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [],
            output: crate::MoleculeBuilder,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn() -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "MoleculeBuilder.from_parts",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::from_parts,
            python: "from_parts",
            javascript: "fromParts",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [
                { name: topology, type: crate::TopologyBlock, default: required },
                { name: coordinates, type: crate::CoordinateBlock, default: required },
                { name: properties, type: crate::MoleculeProperties, default: required },
            ],
            output: crate::MoleculeBuilder,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::TopologyBlock, crate::CoordinateBlock, crate::MoleculeProperties)
                -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "MoleculeBuilder.build",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::build,
            python: "build",
            javascript: "build",
            feature: "runtime", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: none,
            signature: fn(crate::MoleculeBuilder) -> Result<crate::Molecule, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.add_atom",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_atom,
            python: "add_atom",
            javascript: "addAtom",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: atom, type: crate::AtomSpec, default: required }],
            output: crate::AtomId,
            error: none,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::AtomSpec) -> crate::AtomId,
        },
        {
            semantic_id: "MoleculeBuilder.add_bond",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_bond,
            python: "add_bond",
            javascript: "addBond",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: bond, type: crate::BondSpec, default: required }],
            output: crate::BondId,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::BondSpec)
                -> Result<crate::BondId, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.set_atom_formal_charge",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::set_atom_formal_charge,
            python: "set_atom_formal_charge",
            javascript: "setAtomFormalCharge",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [
                { name: atom_id, type: crate::AtomId, default: required },
                { name: formal_charge, type: i8, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::AtomId, i8)
                -> Result<(), crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.set_bond_order",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::set_bond_order,
            python: "set_bond_order",
            javascript: "setBondOrder",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [
                { name: bond_id, type: crate::BondId, default: required },
                { name: order, type: crate::BondOrder, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::BondId, crate::BondOrder)
                -> Result<(), crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.remove_bond_between_atoms",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::remove_bond_between_atoms,
            python: "remove_bond_between_atoms",
            javascript: "removeBondBetweenAtoms",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [
                { name: begin_atom, type: crate::AtomId, default: required },
                { name: end_atom, type: crate::AtomId, default: required },
            ],
            output: bool,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::AtomId, crate::AtomId)
                -> Result<bool, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.degree",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::degree,
            python: "degree",
            javascript: "degree",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: atom_id, type: crate::AtomId, default: required }],
            output: usize,
            error: crate::OperationError,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder, crate::AtomId)
                -> Result<usize, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.neighbor_bonds",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::neighbor_bonds,
            python: "neighbor_bonds",
            javascript: "neighborBonds",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: atom_id, type: crate::AtomId, default: required }],
            output: Vec<crate::BondId>,
            error: crate::OperationError,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder, crate::AtomId)
                -> Result<Vec<crate::BondId>, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.bond_between_atoms",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::bond_between_atoms,
            python: "bond_between_atoms",
            javascript: "bondBetweenAtoms",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [
                { name: begin_atom, type: crate::AtomId, default: required },
                { name: end_atom, type: crate::AtomId, default: required },
            ],
            output: Option<crate::BondId>,
            error: crate::OperationError,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder, crate::AtomId, crate::AtomId)
                -> Result<Option<crate::BondId>, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.atoms",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::atoms,
            python: "atoms",
            javascript: "atoms",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::Atom],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &[crate::Atom],
        },
        {
            semantic_id: "MoleculeBuilder.bonds",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::bonds,
            python: "bonds",
            javascript: "bonds",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::Bond],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &[crate::Bond],
        },
        {
            semantic_id: "MoleculeBuilder.substance_groups",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::substance_groups,
            python: "substance_groups",
            javascript: "substanceGroups",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::SubstanceGroup],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &[crate::SubstanceGroup],
        },
        {
            semantic_id: "MoleculeBuilder.stereo_groups",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::stereo_groups,
            python: "stereo_groups",
            javascript: "stereoGroups",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::StereoGroup],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &[crate::StereoGroup],
        },
        {
            semantic_id: "MoleculeBuilder.coordinates",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::coordinates,
            python: "coordinates",
            javascript: "coordinates",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::CoordinateBlock,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &crate::CoordinateBlock,
        },
        {
            semantic_id: "MoleculeBuilder.properties",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::properties,
            python: "properties",
            javascript: "properties",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::MoleculeProperties,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::MoleculeBuilder) -> &crate::MoleculeProperties,
        },
        {
            semantic_id: "MoleculeBuilder.set_2d_coordinates",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::set_2d_coordinates,
            python: "set_2d_coordinates",
            javascript: "set2dCoordinates",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: coordinates, type: Vec<[f64; 2]>, default: required }],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, Vec<[f64; 2]>)
                -> Result<(), crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.add_2d_conformer",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_2d_conformer,
            python: "add_2d_conformer",
            javascript: "add2dConformer",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: coordinates, type: Vec<[f64; 2]>, default: required }],
            output: usize,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, Vec<[f64; 2]>)
                -> Result<usize, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.add_3d_conformer",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_3d_conformer,
            python: "add_3d_conformer",
            javascript: "add3dConformer",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: coordinates, type: Vec<[f64; 3]>, default: required }],
            output: usize,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, Vec<[f64; 3]>)
                -> Result<usize, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.add_substance_group",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_substance_group,
            python: "add_substance_group",
            javascript: "addSubstanceGroup",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: group, type: crate::SubstanceGroup, default: required }],
            output: crate::SubstanceGroupId,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::SubstanceGroup)
                -> Result<crate::SubstanceGroupId, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.add_stereo_group",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::add_stereo_group,
            python: "add_stereo_group",
            javascript: "addStereoGroup",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [{ name: group, type: crate::StereoGroup, default: required }],
            output: usize,
            error: crate::OperationError,
            state: in_place,
            operation: none,
            signature: fn(&mut crate::MoleculeBuilder, crate::StereoGroup)
                -> Result<usize, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.with_name",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::with_name,
            python: "with_name",
            javascript: "withName",
            feature: "runtime", status: experimental,
            kind: instance, receiver: owned,
            parameters: [{ name: name, type: String, default: required }],
            output: crate::MoleculeBuilder,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::MoleculeBuilder, String) -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "MoleculeBuilder.with_property",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::with_property,
            python: "with_property",
            javascript: "withProperty",
            feature: "runtime", status: experimental,
            kind: instance, receiver: owned,
            parameters: [
                { name: key, type: String, default: required },
                { name: value, type: String, default: required },
            ],
            output: crate::MoleculeBuilder,
            error: crate::OperationError,
            state: value_returning,
            operation: none,
            signature: fn(crate::MoleculeBuilder, String, String)
                -> Result<crate::MoleculeBuilder, crate::OperationError>,
        },
        {
            semantic_id: "MoleculeBuilder.with_sdf_data_field",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::with_sdf_data_field,
            python: "with_sdf_data_field",
            javascript: "withSdfDataField",
            feature: "runtime", status: experimental,
            kind: instance, receiver: owned,
            parameters: [
                { name: key, type: String, default: required },
                { name: value, type: String, default: required },
            ],
            output: crate::MoleculeBuilder,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::MoleculeBuilder, String, String) -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "MoleculeBuilder.with_properties",
            item: callable,
            owner: type_,
            rust: crate::MoleculeBuilder::with_properties,
            python: "with_properties",
            javascript: "withProperties",
            feature: "runtime", status: experimental,
            kind: instance, receiver: owned,
            parameters: [{ name: properties, type: crate::MoleculeProperties, default: required }],
            output: crate::MoleculeBuilder,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::MoleculeBuilder, crate::MoleculeProperties)
                -> crate::MoleculeBuilder,
        },
        {
            semantic_id: "types.SubstanceGroupId",
            item: type,
            owner: type_,
            rust: crate::SubstanceGroupId,
            python: "SubstanceGroupId",
            javascript: "SubstanceGroupId",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.SubstanceGroupKind",
            item: type,
            owner: type_,
            rust: crate::SubstanceGroupKind,
            python: "SubstanceGroupKind",
            javascript: "SubstanceGroupKind",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.SubstanceGroup",
            item: type,
            owner: type_,
            rust: crate::SubstanceGroup,
            python: "SubstanceGroup",
            javascript: "SubstanceGroup",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.SGroupBracket",
            item: type,
            owner: type_,
            rust: crate::SGroupBracket,
            python: "SGroupBracket",
            javascript: "SGroupBracket",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.SGroupCState",
            item: type,
            owner: type_,
            rust: crate::SGroupCState,
            python: "SGroupCState",
            javascript: "SGroupCState",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.SGroupDisplay",
            item: type,
            owner: type_,
            rust: crate::SGroupDisplay,
            python: "SGroupDisplay",
            javascript: "SGroupDisplay",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "SubstanceGroupId.new",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroupId::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [{ name: index, type: usize, default: required }],
            output: crate::SubstanceGroupId,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(usize) -> crate::SubstanceGroupId,
        },
        {
            semantic_id: "SubstanceGroup.new",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [
                { name: id, type: crate::SubstanceGroupId, default: required },
                { name: kind, type: crate::SubstanceGroupKind, default: required },
            ],
            output: crate::SubstanceGroup,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::SubstanceGroupId, crate::SubstanceGroupKind)
                -> crate::SubstanceGroup,
        },
        {
            semantic_id: "SubstanceGroupId.index",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroupId::index,
            python: "index",
            javascript: "index",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SubstanceGroupId) -> usize,
        },
        {
            semantic_id: "SGroupBracket.new",
            item: callable,
            owner: type_,
            rust: crate::SGroupBracket::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [{ name: points, type: [[f64; 3]; 3], default: required }],
            output: crate::SGroupBracket,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn([[f64; 3]; 3]) -> crate::SGroupBracket,
        },
        {
            semantic_id: "SGroupCState.new",
            item: callable,
            owner: type_,
            rust: crate::SGroupCState::new,
            python: "new",
            javascript: "new",
            feature: "runtime", status: experimental,
            kind: static_,
            parameters: [
                { name: bond, type: crate::BondId, default: required },
                { name: vector, type: [f64; 3], default: required },
            ],
            output: crate::SGroupCState,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::BondId, [f64; 3]) -> crate::SGroupCState,
        },
        {
            semantic_id: "SGroupBracket.points",
            item: callable,
            owner: type_,
            rust: crate::SGroupBracket::points,
            python: "points",
            javascript: "points",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[[f64; 3]; 3],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SGroupBracket) -> &'a [[f64; 3]; 3],
        },
        {
            semantic_id: "SGroupCState.bond",
            item: callable,
            owner: type_,
            rust: crate::SGroupCState::bond,
            python: "bond",
            javascript: "bond",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::BondId,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SGroupCState) -> crate::BondId,
        },
        {
            semantic_id: "SGroupCState.vector",
            item: callable,
            owner: type_,
            rust: crate::SGroupCState::vector,
            python: "vector",
            javascript: "vector",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[f64; 3],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SGroupCState) -> &'a [f64; 3],
        },
        {
            semantic_id: "SGroupDisplay.brackets",
            item: callable,
            owner: type_,
            rust: crate::SGroupDisplay::brackets,
            python: "brackets",
            javascript: "brackets",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::SGroupBracket],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SGroupDisplay) -> &'a [crate::SGroupBracket],
        },
        {
            semantic_id: "SubstanceGroup.id",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::id,
            python: "id",
            javascript: "id",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::SubstanceGroupId,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SubstanceGroup) -> crate::SubstanceGroupId,
        },
        {
            semantic_id: "SubstanceGroup.kind",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::kind,
            python: "kind",
            javascript: "kind",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::SubstanceGroupKind,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SubstanceGroup) -> &'a crate::SubstanceGroupKind,
        },
        {
            semantic_id: "SubstanceGroup.display",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::display,
            python: "display",
            javascript: "display",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: Option<&crate::SGroupDisplay>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SubstanceGroup) -> Option<&'a crate::SGroupDisplay>,
        },
        {
            semantic_id: "SubstanceGroup.cstates",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::cstates,
            python: "cstates",
            javascript: "cstates",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::SGroupCState],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SubstanceGroup) -> &'a [crate::SGroupCState],
        },
        {
            semantic_id: "SubstanceGroup.head_crossing_bonds",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::head_crossing_bonds,
            python: "head_crossing_bonds",
            javascript: "headCrossingBonds",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::BondId],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SubstanceGroup) -> &'a [crate::BondId],
        },
        {
            semantic_id: "SubstanceGroup.crossing_bond_correspondence",
            item: callable,
            owner: type_,
            rust: crate::SubstanceGroup::crossing_bond_correspondence,
            python: "crossing_bond_correspondence",
            javascript: "crossingBondCorrespondence",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::BondId],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::SubstanceGroup) -> &'a [crate::BondId],
        },
        {
            semantic_id: "types.OperationError",
            item: type,
            owner: type_,
            rust: crate::OperationError,
            python: "OperationError",
            javascript: "OperationError",
            feature: "runtime", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffParameterQueryError",
            item: type,
            owner: type_,
            rust: crate::UffParameterQueryError,
            python: "UffParameterQueryError",
            javascript: "UffParameterQueryError",
            feature: "cap-forcefields",
            status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffParameterError",
            item: type,
            owner: type_,
            rust: crate::UffParameterError,
            python: "UffParameterError",
            javascript: "UffParameterError",
            feature: "cap-forcefields",
            status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffParameterErrorKind",
            item: type,
            owner: type_,
            rust: crate::UffParameterErrorKind,
            python: "UffParameterErrorKind",
            javascript: "UffParameterErrorKind",
            feature: "cap-forcefields",
            status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "UffParameterError.kind",
            item: callable,
            owner: type_,
            rust: crate::UffParameterError::kind,
            python: "kind",
            javascript: "kind",
            feature: "cap-forcefields",
            status: experimental,
            kind: instance,
            parameters: [],
            output: crate::UffParameterErrorKind,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::UffParameterError) -> crate::UffParameterErrorKind,
        },
        // Canonical detached atom/bond property values. This registry exposes
        // the four modeled variants through their value kind and exact typed
        // borrows only. Source string projection remains a private core/domain
        // conversion boundary and is deliberately not registered here.
        {semantic_id:"types.MoleculeProperties",item:type,owner:type_,rust:crate::MoleculeProperties,python:"MoleculeProperties",javascript:"MoleculeProperties",feature:"metadata",status:experimental,role:value,},
        {semantic_id:"types.SdfPropertyList",item:type,owner:type_,rust:crate::SdfPropertyList,python:"SdfPropertyList",javascript:"SdfPropertyList",feature:"metadata",status:experimental,role:value,},
        {semantic_id:"types.SdfPropertyListTarget",item:type,owner:type_,rust:crate::SdfPropertyListTarget,python:"SdfPropertyListTarget",javascript:"SdfPropertyListTarget",feature:"metadata",status:experimental,role:value,},
        {semantic_id:"MoleculeProperties.name",item:callable,owner:type_,rust:crate::MoleculeProperties::name,python:"name",javascript:"name",feature:"metadata",status:experimental,kind:instance,parameters:[],output:Option<&'a str>,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MoleculeProperties)->Option<&'a str>,},
        {semantic_id:"MoleculeProperties.sdf_data_fields",item:callable,owner:type_,rust:crate::MoleculeProperties::sdf_data_fields,python:"sdf_data_fields",javascript:"sdfDataFields",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a [(String,String)],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MoleculeProperties)->&'a [(String,String)],},
        {semantic_id:"MoleculeProperties.sdf_property_lists",item:callable,owner:type_,rust:crate::MoleculeProperties::sdf_property_lists,python:"sdf_property_lists",javascript:"sdfPropertyLists",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a [crate::SdfPropertyList],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MoleculeProperties)->&'a [crate::SdfPropertyList],},
        {semantic_id:"MoleculeProperties.props",item:callable,owner:type_,rust:crate::MoleculeProperties::props,python:"props",javascript:"props",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a std::collections::BTreeMap<String,String>,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MoleculeProperties)->&'a std::collections::BTreeMap<String,String>,},
        {semantic_id:"MoleculeProperties.prop",item:callable,owner:type_,rust:crate::MoleculeProperties::prop,python:"prop",javascript:"prop",feature:"metadata",status:experimental,kind:instance,parameters:[{name:key,type:&str,default:required}],output:Option<&'a str>,error:none,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::MoleculeProperties,&'b str)->Option<&'a str>,},
        {semantic_id:"MoleculeProperties.is_prop_computed",item:callable,owner:type_,rust:crate::MoleculeProperties::is_prop_computed,python:"is_prop_computed",javascript:"isPropComputed",feature:"metadata",status:experimental,kind:instance,parameters:[{name:key,type:&str,default:required}],output:bool,error:none,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::MoleculeProperties,&'b str)->bool,},
        {semantic_id:"MoleculeProperties.computed_prop_names",item:callable,owner:type_,rust:crate::MoleculeProperties::computed_prop_names,python:"computed_prop_names",javascript:"computedPropNames",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a std::collections::BTreeSet<String>,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::MoleculeProperties)->&'a std::collections::BTreeSet<String>,},
        {semantic_id:"SdfPropertyList.target",item:callable,owner:type_,rust:crate::SdfPropertyList::target,python:"target",javascript:"target",feature:"metadata",status:experimental,kind:instance,parameters:[],output:crate::SdfPropertyListTarget,error:none,state:read_only,operation:none,signature:fn(&crate::SdfPropertyList)->crate::SdfPropertyListTarget,},
        {semantic_id:"SdfPropertyList.name",item:callable,owner:type_,rust:crate::SdfPropertyList::name,python:"name",javascript:"name",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a str,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::SdfPropertyList)->&'a str,},
        {semantic_id:"SdfPropertyList.values",item:callable,owner:type_,rust:crate::SdfPropertyList::values,python:"values",javascript:"values",feature:"metadata",status:experimental,kind:instance,parameters:[],output:&'a [Option<crate::PropertyValue>],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::SdfPropertyList)->&'a [Option<crate::PropertyValue>],},
        {
            semantic_id: "types.PropertyValue",
            item: type,
            owner: type_,
            rust: crate::PropertyValue,
            python: "PropertyValue",
            javascript: "PropertyValue",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.PropertyValueKind",
            item: type,
            owner: type_,
            rust: crate::PropertyValueKind,
            python: "PropertyValueKind",
            javascript: "PropertyValueKind",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.PropertyValueError",
            item: type,
            owner: type_,
            rust: crate::PropertyValueError,
            python: "PropertyValueError",
            javascript: "PropertyValueError",
            feature: "runtime", status: experimental,
            role: error,
        },
        {
            semantic_id: "PropertyValueError.expected",
            item: callable,
            owner: type_,
            rust: crate::PropertyValueError::expected,
            python: "expected",
            javascript: "expected",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::PropertyValueKind,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValueError) -> crate::PropertyValueKind,
        },
        {
            semantic_id: "PropertyValueError.actual",
            item: callable,
            owner: type_,
            rust: crate::PropertyValueError::actual,
            python: "actual",
            javascript: "actual",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::PropertyValueKind,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValueError) -> crate::PropertyValueKind,
        },
        {
            semantic_id: "PropertyValue.kind",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::kind,
            python: "kind",
            javascript: "kind",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::PropertyValueKind,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValue) -> crate::PropertyValueKind,
        },
        {
            semantic_id: "PropertyValue.as_string",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_string,
            python: "as_string",
            javascript: "asString",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &str,
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::PropertyValue,
            ) -> Result<&'a str, crate::PropertyValueError>,
        },
        {
            semantic_id: "PropertyValue.as_int",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_int,
            python: "as_int",
            javascript: "asInt",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: i32,
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValue) -> Result<i32, crate::PropertyValueError>,
        },
        {
            semantic_id: "PropertyValue.as_uint",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_uint,
            python: "as_uint",
            javascript: "asUint",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValue) -> Result<u32, crate::PropertyValueError>,
        },
        {
            semantic_id: "PropertyValue.as_int_vector",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_int_vector,
            python: "as_int_vector",
            javascript: "asIntVector",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[i32],
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::PropertyValue,
            ) -> Result<&'a [i32], crate::PropertyValueError>,
        },
        {
            semantic_id: "PropertyValue.as_double",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_double,
            python: "as_double",
            javascript: "asDouble",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: f64,
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValue) -> Result<f64, crate::PropertyValueError>,
        },
        {
            semantic_id: "PropertyValue.as_bool",
            item: callable,
            owner: type_,
            rust: crate::PropertyValue::as_bool,
            python: "as_bool",
            javascript: "asBool",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: bool,
            error: crate::PropertyValueError,
            state: read_only,
            operation: none,
            signature: fn(&crate::PropertyValue) -> Result<bool, crate::PropertyValueError>,
        },
        // The canonical typed template ATTCHORD state is owned by model and
        // re-exported by cosmolkit. These rows expose only immutable values
        // and accessors; no live-Molecule mutation surface is introduced.
        {
            semantic_id: "types.TemplateAttachment",
            item: type,
            owner: type_,
            rust: crate::TemplateAttachment,
            python: "TemplateAttachment",
            javascript: "TemplateAttachment",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.TemplateAttachmentOrder",
            item: type,
            owner: type_,
            rust: crate::TemplateAttachmentOrder,
            python: "TemplateAttachmentOrder",
            javascript: "TemplateAttachmentOrder",
            feature: "runtime", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.TemplateAttachmentOrderError",
            item: type,
            owner: type_,
            rust: crate::TemplateAttachmentOrderError,
            python: "TemplateAttachmentOrderError",
            javascript: "TemplateAttachmentOrderError",
            feature: "runtime", status: experimental,
            role: error,
        },
        {
            semantic_id: "TemplateAttachment.target",
            item: callable,
            owner: type_,
            rust: crate::TemplateAttachment::target,
            python: "target",
            javascript: "target",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::AtomId,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::TemplateAttachment) -> crate::AtomId,
        },
        {
            semantic_id: "TemplateAttachment.label",
            item: callable,
            owner: type_,
            rust: crate::TemplateAttachment::label,
            python: "label",
            javascript: "label",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &str,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::TemplateAttachment) -> &'a str,
        },
        {
            semantic_id: "TemplateAttachmentOrder.entries",
            item: callable,
            owner: type_,
            rust: crate::TemplateAttachmentOrder::entries,
            python: "entries",
            javascript: "entries",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::TemplateAttachment],
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::TemplateAttachmentOrder,
            ) -> &'a [crate::TemplateAttachment],
        },
        {
            semantic_id: "Atom.template_attachment_order",
            item: callable,
            owner: type_,
            rust: crate::Atom::template_attachment_order,
            python: "template_attachment_order",
            javascript: "templateAttachmentOrder",
            feature: "runtime", status: experimental,
            kind: instance,
            parameters: [],
            output: Option<&crate::TemplateAttachmentOrder>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(
                &'a crate::Atom,
            ) -> Option<&'a crate::TemplateAttachmentOrder>,
        },
        {
            semantic_id: "types.FunctionStatus",
            item: type,
            owner: type_,
            rust: crate::FunctionStatus,
            python: "FunctionStatus",
            javascript: "FunctionStatus",
            feature: "metadata", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.ParityPolicy",
            item: type,
            owner: type_,
            rust: crate::ParityPolicy,
            python: "ParityPolicy",
            javascript: "ParityPolicy",
            feature: "metadata", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.FeatureSpec",
            item: type,
            owner: type_,
            rust: crate::FeatureSpec,
            python: "FeatureSpec",
            javascript: "FeatureSpec",
            feature: "metadata", status: experimental,
            role: value,
        },
        {
            semantic_id: "types.MoleculeOpSpec",
            item: type,
            owner: type_,
            rust: crate::MoleculeOpSpec,
            python: "MoleculeOpSpec",
            javascript: "MoleculeOpSpec",
            feature: "metadata", status: experimental,
            role: result,
        },
        {
            semantic_id: "types.SupportMatrixEntry",
            item: type,
            owner: type_,
            rust: crate::SupportMatrixEntry,
            python: "SupportMatrixEntry",
            javascript: "SupportMatrixEntry",
            feature: "metadata", status: experimental,
            role: result,
        },
        {
            semantic_id: "types.OperationInvariantEntry",
            item: type,
            owner: type_,
            rust: crate::OperationInvariantEntry,
            python: "OperationInvariantEntry",
            javascript: "OperationInvariantEntry",
            feature: "metadata", status: experimental,
            role: result,
        },
        {
            semantic_id: "types.ParityMatrixEntry",
            item: type,
            owner: type_,
            rust: crate::ParityMatrixEntry,
            python: "ParityMatrixEntry",
            javascript: "ParityMatrixEntry",
            feature: "metadata", status: experimental,
            role: result,
        },
        {
            semantic_id: "types.FeatureSpecIter",
            item: type,
            owner: type_,
            rust: crate::FeatureSpecIter,
            python: "FeatureSpecIter",
            javascript: "FeatureSpecIter",
            feature: "metadata", status: experimental,
            role: result,
        },
        {
            semantic_id: "module.feature_specs",
            item: callable,
            owner: module,
            rust: crate::feature_specs,
            python: "feature_specs",
            javascript: "featureSpecs",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: crate::FeatureSpecIter,
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> crate::FeatureSpecIter,
        },
        {
            semantic_id: "module.feature_spec",
            item: callable,
            owner: module,
            rust: crate::feature_spec,
            python: "feature_spec",
            javascript: "featureSpec",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [{ name: name, type: &str, default: required }],
            output: Option<&'static crate::FeatureSpec>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> Option<&'static crate::FeatureSpec>,
        },
        {
            semantic_id: "module.operation_specs",
            item: callable,
            owner: module,
            rust: crate::operation_specs,
            python: "operation_specs",
            javascript: "operationSpecs",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: &'static [&'static crate::MoleculeOpSpec],
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> &'static [&'static crate::MoleculeOpSpec],
        },
        {
            semantic_id: "module.operation_spec",
            item: callable,
            owner: module,
            rust: crate::operation_spec,
            python: "operation_spec",
            javascript: "operationSpec",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [{ name: method, type: &str, default: required }],
            output: Option<&'static crate::MoleculeOpSpec>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> Option<&'static crate::MoleculeOpSpec>,
        },
        {
            semantic_id: "module.support_matrix",
            item: callable,
            owner: module,
            rust: crate::support_matrix,
            python: "support_matrix",
            javascript: "supportMatrix",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: &'static [crate::SupportMatrixEntry],
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> &'static [crate::SupportMatrixEntry],
        },
        {
            semantic_id: "module.operation_invariant_matrix",
            item: callable,
            owner: module,
            rust: crate::operation_invariant_matrix,
            python: "operation_invariant_matrix",
            javascript: "operationInvariantMatrix",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: &'static [crate::OperationInvariantEntry],
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> &'static [crate::OperationInvariantEntry],
        },
        {
            semantic_id: "module.parity_matrix",
            item: callable,
            owner: module,
            rust: crate::parity_matrix,
            python: "parity_matrix",
            javascript: "parityMatrix",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: &'static [crate::ParityMatrixEntry],
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> &'static [crate::ParityMatrixEntry],
        },
        {
            semantic_id: "module.operation_invariant",
            item: callable,
            owner: module,
            rust: crate::operation_invariant,
            python: "operation_invariant",
            javascript: "operationInvariant",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [{ name: method, type: &str, default: required }],
            output: Option<&'static crate::OperationInvariantEntry>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> Option<&'static crate::OperationInvariantEntry>,
        },
        {
            semantic_id: "module.operation_parity",
            item: callable,
            owner: module,
            rust: crate::operation_parity,
            python: "operation_parity",
            javascript: "operationParity",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [{ name: method, type: &str, default: required }],
            output: Option<&'static crate::ParityMatrixEntry>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> Option<&'static crate::ParityMatrixEntry>,
        },
        {
            semantic_id: "module.version",
            item: callable,
            owner: module,
            rust: crate::version,
            python: "version",
            javascript: "version",
            feature: "metadata", status: experimental,
            kind: module,
            parameters: [],
            output: &'static str,
            error: none,
            state: read_only,
            operation: none,
            signature: fn() -> &'static str,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "types.DenseMatrix",
            item: type,
            owner: type_,
            rust: crate::DenseMatrix,
            python: "DenseMatrix",
            javascript: "DenseMatrix",
            feature: "cap-matrices", status: experimental,
            role: result,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "types.DistanceMatrixParams",
            item: type,
            owner: type_,
            rust: crate::DistanceMatrixParams,
            python: "DistanceMatrixParams",
            javascript: "DistanceMatrixParams",
            feature: "cap-matrices", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "types.DistanceMatrix3dParams",
            item: type,
            owner: type_,
            rust: crate::DistanceMatrix3dParams,
            python: "DistanceMatrix3dParams",
            javascript: "DistanceMatrix3dParams",
            feature: "cap-matrices", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "types.MatrixError",
            item: type,
            owner: type_,
            rust: crate::MatrixError,
            python: "MatrixError",
            javascript: "MatrixError",
            feature: "cap-matrices", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "Molecule.distance_matrix",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::distance_matrix,
            python: "distance_matrix",
            javascript: "distanceMatrix",
            feature: "cap-matrices", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::DenseMatrix,
            error: crate::MatrixError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<crate::DenseMatrix, crate::MatrixError>,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "Molecule.distance_matrix_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::distance_matrix_with_params,
            python: "distance_matrix_with_params",
            javascript: "distanceMatrixWithParams",
            feature: "cap-matrices", status: experimental,
            kind: instance,
            parameters: [
                { name: params, type: &crate::DistanceMatrixParams, default: required },
            ],
            output: crate::DenseMatrix,
            error: crate::MatrixError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::DistanceMatrixParams,
            ) -> Result<crate::DenseMatrix, crate::MatrixError>,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "Molecule.distance_matrix_3d",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::distance_matrix_3d,
            python: "distance_matrix_3d",
            javascript: "distanceMatrix3d",
            feature: "cap-matrices", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::DenseMatrix,
            error: crate::MatrixError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<crate::DenseMatrix, crate::MatrixError>,
        },
        #[cfg(feature = "cap-matrices")]
        {
            semantic_id: "Molecule.distance_matrix_3d_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::distance_matrix_3d_with_params,
            python: "distance_matrix_3d_with_params",
            javascript: "distanceMatrix3dWithParams",
            feature: "cap-matrices", status: experimental,
            kind: instance,
            parameters: [
                { name: params, type: &crate::DistanceMatrix3dParams, default: required },
            ],
            output: crate::DenseMatrix,
            error: crate::MatrixError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::DistanceMatrix3dParams,
            ) -> Result<crate::DenseMatrix, crate::MatrixError>,
        },
        // DRAW-layout delegates through the generated coordinate operation to
        // the unique detached depict owner. The short method
        // uses Coordinate2DParams::default(): no constraints, no canonical
        // orientation, replace current 2D rows, zero samples/flips/seed,
        // no degree-four permutations, forceRDKit=false, no ring templates.
        // The configured method carries explicit atom-indexed XY constraints;
        // width/height are renderer options and do not belong here. Both
        // methods are value-returning coordinate operations: read topology,
        // write only the coordinate block and its derived-state bookkeeping,
        // preserve independent 3D rows and source-object identity, with no
        // topology edit or mapping. The generated operation declaration owns
        // their shared function status as well as this execution contract.
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.Coordinate2DParams",
            item: type,
            owner: type_,
            rust: crate::Coordinate2DParams,
            python: "Coordinate2DParams",
            javascript: "Coordinate2DParams",
            feature: "cap-depict", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.Coordinate2DError",
            item: type,
            owner: type_,
            rust: crate::Coordinate2DError,
            python: "Coordinate2DError",
            javascript: "Coordinate2DError",
            feature: "cap-depict", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.Coordinate2DTemplateError",
            item: type,
            owner: type_,
            rust: crate::Coordinate2DTemplateError,
            python: "Coordinate2DTemplateError",
            javascript: "Coordinate2DTemplateError",
            feature: "cap-depict", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.Coordinate2DLayoutError",
            item: type,
            owner: type_,
            rust: crate::Coordinate2DLayoutError,
            python: "Coordinate2DLayoutError",
            javascript: "Coordinate2DLayoutError",
            feature: "cap-depict", status: experimental,
            role: error,
        },
        {
            semantic_id: "Molecule.has_2d_coordinates", item: callable, owner: molecule,
            rust: crate::Molecule::has_2d_coordinates, python: "has_2d_coordinates", javascript: "has2dCoordinates",
            feature: "runtime", status: experimental, kind: instance,
            parameters: [], output: bool, error: none, state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> bool,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.with_2d_coordinates",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_2d_coordinates,
            python: "with_2d_coordinates",
            javascript: "with2dCoordinates",
            feature: "cap-depict",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_2d_coordinates",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.with_2d_coordinates_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_2d_coordinates_with_params,
            python: "with_2d_coordinates_with_params",
            javascript: "with2dCoordinatesWithParams",
            feature: "cap-depict",
            kind: instance,
            parameters: [
                { name: params, type: &crate::Coordinate2DParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_2d_coordinates_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::Coordinate2DParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.DrawingError", item: type, owner: type_,
            rust: crate::DrawingError, python: "DrawingError", javascript: "DrawingError",
            feature: "cap-depict", status: experimental, role: error,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.to_svg", item: callable, owner: molecule,
            rust: crate::Molecule::to_svg, python: "to_svg", javascript: "toSvg",
            feature: "cap-depict",
            status: parity_with_differences(
                "RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 MolDraw2DSVG",
                "ROOT-SVG-CANONICAL-METADATA-20261005: public SVG declares ck=https://kit.cosmol.org/ instead of the pinned source renderer identity; all other drawing bytes retain their source comparison."
            ), kind: instance,
            parameters: [
                { name: width, type: u32, default: required },
                { name: height, type: u32, default: required },
            ],
            output: String, error: crate::DrawingError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, u32, u32) -> Result<String, crate::DrawingError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.to_png", item: callable, owner: molecule,
            rust: crate::Molecule::to_png, python: "to_png", javascript: "toPng",
            feature: "cap-depict", status: experimental, kind: instance,
            parameters: [
                { name: width, type: u32, default: required },
                { name: height, type: u32, default: required },
            ],
            output: Vec<u8>, error: crate::DrawingError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, u32, u32) -> Result<Vec<u8>, crate::DrawingError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.compute_2d_coordinates_", item: callable, owner: molecule,
            rust: crate::Molecule::compute_2d_coordinates_, python: "compute_2d_coordinates_", javascript: "compute2dCoordinates",
            feature: "cap-depict", kind: instance,
            parameters: [], output: (), error: crate::OperationError,
            state: in_place, operation: "compute_2d_coordinates_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.compute_2d_coordinates_with_params_", item: callable, owner: molecule,
            rust: crate::Molecule::compute_2d_coordinates_with_params_, python: "compute_2d_coordinates_with_params_", javascript: "compute2dCoordinatesWithParams",
            feature: "cap-depict", kind: instance,
            parameters: [{ name: params, type: &crate::Coordinate2DParams, default: required }], output: (), error: crate::OperationError,
            state: in_place, operation: "compute_2d_coordinates_with_params_",
            signature: for<'a, 'b> fn(&'a mut crate::Molecule, &'b crate::Coordinate2DParams) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "types.DrawingWriteError", item: type, owner: type_,
            rust: crate::DrawingWriteError, python: "DrawingWriteError", javascript: "DrawingWriteError",
            feature: "cap-depict", status: experimental, role: error,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.write_svg", item: callable, owner: molecule,
            rust: crate::Molecule::write_svg, python: "write_svg", javascript: "writeSvg",
            feature: "cap-depict",
            status: parity_with_differences(
                "RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 MolDraw2DSVG",
                "ROOT-SVG-CANONICAL-METADATA-20261005: public SVG declares ck=https://kit.cosmol.org/ instead of the pinned source renderer identity; all other drawing bytes retain their source comparison."
            ), kind: instance,
            parameters: [
                { name: path, type: &std::path::Path, default: required },
                { name: width, type: u32, default: required },
                { name: height, type: u32, default: required },
            ],
            output: (), error: crate::DrawingWriteError, state: read_only, operation: none,
            signature: fn(&crate::Molecule, &std::path::Path, u32, u32) -> Result<(), crate::DrawingWriteError>,
        },
        #[cfg(feature = "cap-depict")]
        {
            semantic_id: "Molecule.write_png", item: callable, owner: molecule,
            rust: crate::Molecule::write_png, python: "write_png", javascript: "writePng",
            feature: "cap-depict", status: experimental, kind: instance,
            parameters: [
                { name: path, type: &std::path::Path, default: required },
                { name: width, type: u32, default: required },
                { name: height, type: u32, default: required },
            ],
            output: (), error: crate::DrawingWriteError, state: read_only, operation: none,
            signature: fn(&crate::Molecule, &std::path::Path, u32, u32) -> Result<(), crate::DrawingWriteError>,
        },
        // CORE-transforms public projection delegates to the single detached
        // owner through the generated coordinate-operation runtime.
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "types.AtomPositionParams",
            item: type,
            owner: type_,
            rust: crate::AtomPositionParams,
            python: "AtomPositionParams",
            javascript: "AtomPositionParams",
            feature: "cap-transforms", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "types.TransformError",
            item: type,
            owner: type_,
            rust: crate::TransformError,
            python: "TransformError",
            javascript: "TransformError",
            feature: "cap-transforms", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "Molecule.with_atom_position",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_atom_position,
            python: "with_atom_position",
            javascript: "withAtomPosition",
            feature: "cap-transforms",
            kind: instance,
            parameters: [
                { name: atom, type: crate::AtomId, default: required },
                { name: position, type: [f64; 3], default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_atom_position",
            signature: fn(
                &crate::Molecule,
                crate::AtomId,
                [f64; 3],
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "Molecule.with_atom_position_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_atom_position_with_params,
            python: "with_atom_position_with_params",
            javascript: "withAtomPositionWithParams",
            feature: "cap-transforms",
            kind: instance,
            parameters: [
                { name: atom, type: crate::AtomId, default: required },
                { name: position, type: [f64; 3], default: required },
                { name: params, type: &crate::AtomPositionParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_atom_position_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                crate::AtomId,
                [f64; 3],
                &'b crate::AtomPositionParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "Molecule.set_atom_position_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::set_atom_position_,
            python: "set_atom_position_",
            javascript: "setAtomPosition",
            feature: "cap-transforms",
            kind: instance,
            parameters: [
                { name: atom, type: crate::AtomId, default: required },
                { name: position, type: [f64; 3], default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "set_atom_position_",
            signature: fn(
                &mut crate::Molecule,
                crate::AtomId,
                [f64; 3],
            ) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-transforms")]
        {
            semantic_id: "Molecule.set_atom_position_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::set_atom_position_with_params_,
            python: "set_atom_position_with_params_",
            javascript: "setAtomPositionWithParams",
            feature: "cap-transforms",
            kind: instance,
            parameters: [
                { name: atom, type: crate::AtomId, default: required },
                { name: position, type: [f64; 3], default: required },
                { name: params, type: &crate::AtomPositionParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "set_atom_position_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                crate::AtomId,
                [f64; 3],
                &'b crate::AtomPositionParams,
            ) -> Result<(), crate::OperationError>,
        },
        // CHEM-sanitize projects the unique detached core owner through one
        // generated mutation operation plus two read-only diagnostics.
        // `SanitizeOperations` owns the exact source flag vocabulary and
        // `SanitizeStage` owns the failed-operation vocabulary; individual
        // constants/variants are not parallel callable registry entries.
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.SanitizeOperations",
            item: type,
            owner: type_,
            rust: crate::SanitizeOperations,
            python: "SanitizeOperations",
            javascript: "SanitizeOperations",
            feature: "cap-sanitize", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.SanitizeStage",
            item: type,
            owner: type_,
            rust: crate::SanitizeStage,
            python: "SanitizeStage",
            javascript: "SanitizeStage",
            feature: "cap-sanitize", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.SanitizeParams",
            item: type,
            owner: type_,
            rust: crate::SanitizeParams,
            python: "SanitizeParams",
            javascript: "SanitizeParams",
            feature: "cap-sanitize", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.SanitizeError",
            item: type,
            owner: type_,
            rust: crate::SanitizeError,
            python: "SanitizeError",
            javascript: "SanitizeError",
            feature: "cap-sanitize", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.ChemistryProblemError",
            item: type,
            owner: type_,
            rust: crate::ChemistryProblemError,
            python: "ChemistryProblemError",
            javascript: "ChemistryProblemError",
            feature: "cap-sanitize", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.ChemistryProblem",
            item: type,
            owner: type_,
            rust: crate::ChemistryProblem,
            python: "ChemistryProblem",
            javascript: "ChemistryProblem",
            feature: "cap-sanitize", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "types.ChemistryProblemReport",
            item: type,
            owner: type_,
            rust: crate::ChemistryProblemReport,
            python: "ChemistryProblemReport",
            javascript: "ChemistryProblemReport",
            feature: "cap-sanitize", status: experimental,
            role: result,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "Molecule.sanitize",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::sanitize,
            python: "sanitize",
            javascript: "sanitize",
            feature: "cap-sanitize",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "sanitize",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "Molecule.sanitize_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::sanitize_with_params,
            python: "sanitize_with_params",
            javascript: "sanitizeWithParams",
            feature: "cap-sanitize",
            kind: instance,
            parameters: [
                { name: params, type: &crate::SanitizeParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "sanitize_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::SanitizeParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "Molecule.detect_chemistry_problems",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::detect_chemistry_problems,
            python: "detect_chemistry_problems",
            javascript: "detectChemistryProblems",
            feature: "cap-sanitize", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::ChemistryProblemReport,
            error: crate::SanitizeError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule)
                -> Result<crate::ChemistryProblemReport, crate::SanitizeError>,
        },
        #[cfg(feature = "cap-sanitize")]
        {
            semantic_id: "Molecule.detect_chemistry_problems_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::detect_chemistry_problems_with_params,
            python: "detect_chemistry_problems_with_params",
            javascript: "detectChemistryProblemsWithParams",
            feature: "cap-sanitize", status: experimental,
            kind: instance,
            parameters: [
                { name: params, type: &crate::SanitizeParams, default: required },
            ],
            output: crate::ChemistryProblemReport,
            error: crate::SanitizeError,
            state: read_only,
            operation: none,
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::SanitizeParams,
            ) -> Result<crate::ChemistryProblemReport, crate::SanitizeError>,
        },
        // CHEM-kekulize is projected through the generated operation runtime;
        // the detached behavior remains uniquely owned by cosmolkit-core.
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "types.KekulizeParams",
            item: type,
            owner: type_,
            rust: crate::KekulizeParams,
            python: "KekulizeParams",
            javascript: "KekulizeParams",
            feature: "cap-kekulize", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "types.KekulizeError",
            item: type,
            owner: type_,
            rust: crate::KekulizeError,
            python: "KekulizeError",
            javascript: "KekulizeError",
            feature: "cap-kekulize", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "Molecule.with_kekulized_bonds",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_kekulized_bonds,
            python: "with_kekulized_bonds",
            javascript: "withKekulizedBonds",
            feature: "cap-kekulize",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_kekulized_bonds",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "Molecule.with_kekulized_bonds_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_kekulized_bonds_with_params,
            python: "with_kekulized_bonds_with_params",
            javascript: "withKekulizedBondsWithParams",
            feature: "cap-kekulize",
            kind: instance,
            parameters: [
                { name: params, type: &crate::KekulizeParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_kekulized_bonds_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::KekulizeParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "Molecule.kekulize_bonds_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::kekulize_bonds_,
            python: "kekulize_bonds_",
            javascript: "kekulizeBonds",
            feature: "cap-kekulize",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "kekulize_bonds_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-kekulize")]
        {
            semantic_id: "Molecule.kekulize_bonds_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::kekulize_bonds_with_params_,
            python: "kekulize_bonds_with_params_",
            javascript: "kekulizeBondsWithParams",
            feature: "cap-kekulize",
            kind: instance,
            parameters: [
                { name: params, type: &crate::KekulizeParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "kekulize_bonds_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::KekulizeParams,
            ) -> Result<(), crate::OperationError>,
        },
        // CHEM-aromaticity is projected through the generated operation runtime;
        // the source-backed behavior remains uniquely owned by cosmolkit-core.
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "types.AromaticityModel",
            item: type,
            owner: type_,
            rust: crate::AromaticityModel,
            python: "AromaticityModel",
            javascript: "AromaticityModel",
            feature: "cap-aromaticity", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "types.AromaticityParams",
            item: type,
            owner: type_,
            rust: crate::AromaticityParams,
            python: "AromaticityParams",
            javascript: "AromaticityParams",
            feature: "cap-aromaticity", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "types.AromaticityError",
            item: type,
            owner: type_,
            rust: crate::AromaticityError,
            python: "AromaticityError",
            javascript: "AromaticityError",
            feature: "cap-aromaticity", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "Molecule.with_assigned_aromaticity",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_aromaticity,
            python: "with_assigned_aromaticity",
            javascript: "withAssignedAromaticity",
            feature: "cap-aromaticity",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_aromaticity",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "Molecule.with_assigned_aromaticity_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_aromaticity_with_params,
            python: "with_assigned_aromaticity_with_params",
            javascript: "withAssignedAromaticityWithParams",
            feature: "cap-aromaticity",
            kind: instance,
            parameters: [
                { name: params, type: &crate::AromaticityParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_aromaticity_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::AromaticityParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "Molecule.assign_aromaticity_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_aromaticity_,
            python: "assign_aromaticity_",
            javascript: "assignAromaticity",
            feature: "cap-aromaticity",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_aromaticity_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-aromaticity")]
        {
            semantic_id: "Molecule.assign_aromaticity_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_aromaticity_with_params_,
            python: "assign_aromaticity_with_params_",
            javascript: "assignAromaticityWithParams",
            feature: "cap-aromaticity",
            kind: instance,
            parameters: [
                { name: params, type: &crate::AromaticityParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_aromaticity_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::AromaticityParams,
            ) -> Result<(), crate::OperationError>,
        },
        // CORE-valence public projection: these frozen names and signatures
        // delegate to the single detached owner through the live runtime.
        // The mutation projections share the single `with_assigned_valence`
        // operation declaration; this registry does not define an implementation.
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "types.ValenceModel",
            item: type,
            owner: type_,
            rust: crate::ValenceModel,
            python: "ValenceModel",
            javascript: "ValenceModel",
            feature: "cap-valence", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "types.ValenceParams",
            item: type,
            owner: type_,
            rust: crate::ValenceParams,
            python: "ValenceParams",
            javascript: "ValenceParams",
            feature: "cap-valence", status: experimental,
            role: parameter,
        },
                {
            semantic_id: "types.ValenceError",
            item: type,
            owner: type_,
            rust: crate::ValenceError,
            python: "ValenceError",
            javascript: "ValenceError",
            feature: "metadata", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "Molecule.with_assigned_valence",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_valence,
            python: "with_assigned_valence",
            javascript: "withAssignedValence",
            feature: "cap-valence",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_valence",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "Molecule.with_assigned_valence_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_valence_with_params,
            python: "with_assigned_valence_with_params",
            javascript: "withAssignedValenceWithParams",
            feature: "cap-valence",
            kind: instance,
            parameters: [{ name: params, type: &crate::ValenceParams, default: required }],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_valence_with_params",
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::ValenceParams)
                -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "Molecule.assign_valence_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_valence_,
            python: "assign_valence_",
            javascript: "assignValence",
            feature: "cap-valence",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_valence_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "Molecule.assign_valence_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_valence_with_params_,
            python: "assign_valence_with_params_",
            javascript: "assignValenceWithParams",
            feature: "cap-valence",
            kind: instance,
            parameters: [{ name: params, type: &crate::ValenceParams, default: required }],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_valence_with_params_",
            signature: for<'a, 'b> fn(&'a mut crate::Molecule, &'b crate::ValenceParams)
                -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-valence")]
        {
            semantic_id: "Molecule.has_valence_violation",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::has_valence_violation,
            python: "has_valence_violation",
            javascript: "hasValenceViolation",
            feature: "cap-valence", status: experimental,
            kind: instance,
            parameters: [{ name: atom_id, type: crate::AtomId, default: required }],
            output: bool,
            error: crate::ValenceError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule, crate::AtomId)
                -> Result<bool, crate::ValenceError>,
        },
        #[cfg(feature = "cap-radicals")]
        {
            semantic_id: "Molecule.with_assigned_radicals",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_radicals,
            python: "with_assigned_radicals",
            javascript: "withAssignedRadicals",
            feature: "cap-radicals",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_radicals",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-radicals")]
        {
            semantic_id: "Molecule.assign_radicals_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_radicals_,
            python: "assign_radicals_",
            javascript: "assignRadicals",
            feature: "cap-radicals",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_radicals_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "types.RingSearchParams",
            item: type,
            owner: type_,
            rust: crate::RingSearchParams,
            python: "RingSearchParams",
            javascript: "RingSearchParams",
            feature: "cap-rings", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.with_assigned_rings",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_rings,
            python: "with_assigned_rings",
            javascript: "withAssignedRings",
            feature: "cap-rings",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_rings",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.assign_rings_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_rings_,
            python: "assign_rings_",
            javascript: "assignRings",
            feature: "cap-rings",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_rings_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.with_assigned_ring_families",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_ring_families,
            python: "with_assigned_ring_families",
            javascript: "withAssignedRingFamilies",
            feature: "cap-rings",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_ring_families",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.with_assigned_ring_families_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_assigned_ring_families_with_params,
            python: "with_assigned_ring_families_with_params",
            javascript: "withAssignedRingFamiliesWithParams",
            feature: "cap-rings",
            kind: instance,
            parameters: [{ name: params, type: &crate::RingSearchParams, default: required }],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_assigned_ring_families_with_params",
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::RingSearchParams)
                -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.assign_ring_families_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_ring_families_,
            python: "assign_ring_families_",
            javascript: "assignRingFamilies",
            feature: "cap-rings",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_ring_families_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-rings")]
        {
            semantic_id: "Molecule.assign_ring_families_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_ring_families_with_params_,
            python: "assign_ring_families_with_params_",
            javascript: "assignRingFamiliesWithParams",
            feature: "cap-rings",
            kind: instance,
            parameters: [{ name: params, type: &crate::RingSearchParams, default: required }],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_ring_families_with_params_",
            signature: for<'a, 'b> fn(&'a mut crate::Molecule, &'b crate::RingSearchParams)
                -> Result<(), crate::OperationError>,
        },
        // CORE-structure_tags exposes one source-backed operation family. The
        // short variants select the source
        // defaults (`conformer_id = -1`, `replace_existing_tags = true`); the
        // `_with_params` variants expose exactly those two fields through one
        // immutable parameter value. The public result is the canonical
        // Molecule/() operation result, never the detached assignment block.
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.StructureTagParams",
            item: type,
            owner: type_,
            rust: crate::StructureTagParams,
            python: "StructureTagParams",
            javascript: "StructureTagParams",
            feature: "cap-stereo", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.StereoError",
            item: type,
            owner: type_,
            rust: crate::StereoError,
            python: "StereoError",
            javascript: "StereoError",
            feature: "cap-stereo", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.with_chiral_tags_from_structure",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_chiral_tags_from_structure,
            python: "with_chiral_tags_from_structure",
            javascript: "withChiralTagsFromStructure",
            feature: "cap-stereo",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_chiral_tags_from_structure",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.with_chiral_tags_from_structure_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_chiral_tags_from_structure_with_params,
            python: "with_chiral_tags_from_structure_with_params",
            javascript: "withChiralTagsFromStructureWithParams",
            feature: "cap-stereo",
            kind: instance,
            parameters: [
                { name: params, type: &crate::StructureTagParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_chiral_tags_from_structure_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::StructureTagParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.assign_chiral_tags_from_structure_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_chiral_tags_from_structure_,
            python: "assign_chiral_tags_from_structure_",
            javascript: "assignChiralTagsFromStructure",
            feature: "cap-stereo",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_chiral_tags_from_structure_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.assign_chiral_tags_from_structure_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_chiral_tags_from_structure_with_params_,
            python: "assign_chiral_tags_from_structure_with_params_",
            javascript: "assignChiralTagsFromStructureWithParams",
            feature: "cap-stereo",
            kind: instance,
            parameters: [
                { name: params, type: &crate::StructureTagParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_chiral_tags_from_structure_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::StructureTagParams,
            ) -> Result<(), crate::OperationError>,
        },
        // CORE-potential is registered before its live wrapper is enabled.
        // These rows freeze the canonical public spelling and typed result
        // surface; the source-backed perception algorithm remains uniquely
        // owned by `cosmolkit-core`.
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoParams",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoParams,
            python: "PotentialStereoParams",
            javascript: "PotentialStereoParams",
            feature: "cap-stereo", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoType",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoType,
            python: "PotentialStereoType",
            javascript: "PotentialStereoType",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoSpecified",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoSpecified,
            python: "PotentialStereoSpecified",
            javascript: "PotentialStereoSpecified",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoDescriptor",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoDescriptor,
            python: "PotentialStereoDescriptor",
            javascript: "PotentialStereoDescriptor",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoCenter",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoCenter,
            python: "PotentialStereoCenter",
            javascript: "PotentialStereoCenter",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoInfo",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoInfo,
            python: "PotentialStereoInfo",
            javascript: "PotentialStereoInfo",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.RingStereoRelation",
            item: type,
            owner: type_,
            rust: crate::RingStereoRelation,
            python: "RingStereoRelation",
            javascript: "RingStereoRelation",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            // The binding surface is the default Molecule specialization.
            // PendingMolecule is private transaction state, never a binding type.
            semantic_id: "types.PotentialStereoResult",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoResult,
            python: "PotentialStereoResult",
            javascript: "PotentialStereoResult",
            feature: "cap-stereo",
            status: experimental,
            role: result,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.PotentialStereoError",
            item: type,
            owner: type_,
            rust: crate::PotentialStereoError,
            python: "PotentialStereoError",
            javascript: "PotentialStereoError",
            feature: "cap-stereo", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.potential_stereo",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::potential_stereo,
            python: "potential_stereo",
            javascript: "potentialStereo",
            feature: "cap-stereo",
            kind: instance,
            parameters: [],
            output: crate::PotentialStereoResult,
            error: crate::OperationError,
            state: read_only,
            operation: "potential_stereo",
            signature: fn(&crate::Molecule)
                -> Result<crate::PotentialStereoResult, crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.potential_stereo_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::potential_stereo_with_params,
            python: "potential_stereo_with_params",
            javascript: "potentialStereoWithParams",
            feature: "cap-stereo",
            kind: instance,
            parameters: [
                { name: params, type: &crate::PotentialStereoParams, default: required },
            ],
            output: crate::PotentialStereoResult,
            error: crate::OperationError,
            state: value_returning,
            operation: "potential_stereo_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::PotentialStereoParams,
            ) -> Result<crate::PotentialStereoResult, crate::OperationError>,
        },
        // ST-cip_labels canonical surface projects the sole detached stereo
        // owner through the registered parent transaction. `CipDescriptor`
        // itself is the complete detached
        // twelve-spelling value vocabulary R/S/r/s/E/Z/e/z/M/P/m/p; in
        // particular RDKit seqTrans/seqCis project to typed LowerE/LowerZ, not
        // invalid stored values. The detached CipLabelAssignment is
        // transaction hand-off state, not a binding result type.
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.CipDescriptor",
            item: type,
            owner: type_,
            rust: crate::CipDescriptor,
            python: "CipDescriptor",
            javascript: "CipDescriptor",
            feature: "cap-stereo", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.CipDescriptorError",
            item: type,
            owner: type_,
            rust: crate::CipDescriptorError,
            python: "CipDescriptorError",
            javascript: "CipDescriptorError",
            feature: "cap-stereo", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.CipLabelOptions",
            item: type,
            owner: type_,
            rust: crate::CipLabelOptions,
            python: "CipLabelOptions",
            javascript: "CipLabelOptions",
            feature: "cap-stereo", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "types.CipLabelerError",
            item: type,
            owner: type_,
            rust: crate::CipLabelerError,
            python: "CipLabelerError",
            javascript: "CipLabelerError",
            feature: "cap-stereo", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.with_cip_labels",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_cip_labels,
            python: "with_cip_labels",
            javascript: "withCipLabels",
            feature: "cap-stereo",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_cip_labels",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.with_cip_labels_with_options",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_cip_labels_with_options,
            python: "with_cip_labels_with_options",
            javascript: "withCipLabelsWithOptions",
            feature: "cap-stereo",
            kind: instance,
            parameters: [
                { name: options, type: &crate::CipLabelOptions, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_cip_labels_with_options",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::CipLabelOptions,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.assign_cip_labels_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_cip_labels_,
            python: "assign_cip_labels_",
            javascript: "assignCipLabels",
            feature: "cap-stereo",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_cip_labels_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-stereo")]
        {
            semantic_id: "Molecule.assign_cip_labels_with_options_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::assign_cip_labels_with_options_,
            python: "assign_cip_labels_with_options_",
            javascript: "assignCipLabelsWithOptions",
            feature: "cap-stereo",
            kind: instance,
            parameters: [
                { name: options, type: &crate::CipLabelOptions, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "assign_cip_labels_with_options_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::CipLabelOptions,
            ) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_0",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_0,
            python: "chi_0", javascript: "chi0",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_1",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_1,
            python: "chi_1", javascript: "chi1",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        // D02 canonical recomputation-only read queries; no operation lifecycle.
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.hall_kier_alpha",
            item: callable, owner: molecule,
            rust: crate::Molecule::hall_kier_alpha,
            python: "hall_kier_alpha", javascript: "hallKierAlpha",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.hall_kier_alpha_with_contributions",
            item: callable, owner: molecule,
            rust: crate::Molecule::hall_kier_alpha_with_contributions,
            python: "hall_kier_alpha_with_contributions", javascript: "hallKierAlphaWithContributions",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: (f64, Vec<f64>), error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<(f64, Vec<f64>), crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.kappa_1",
            item: callable, owner: molecule,
            rust: crate::Molecule::kappa_1,
            python: "kappa_1", javascript: "kappa1",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.kappa_2",
            item: callable, owner: molecule,
            rust: crate::Molecule::kappa_2,
            python: "kappa_2", javascript: "kappa2",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.kappa_3",
            item: callable, owner: molecule,
            rust: crate::Molecule::kappa_3,
            python: "kappa_3", javascript: "kappa3",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.phi",
            item: callable, owner: molecule,
            rust: crate::Molecule::phi,
            python: "phi", javascript: "phi",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.mqns",
            item: callable, owner: molecule,
            rust: crate::Molecule::mqns,
            python: "mqns", javascript: "mqns",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: force, type: bool, default: false }],
            output: Vec<u32>, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool) -> Result<Vec<u32>, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_0_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_0_v,
            python: "chi_0_v", javascript: "chi0V",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_1_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_1_v,
            python: "chi_1_v", javascript: "chi1V",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_2_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_2_v,
            python: "chi_2_v", javascript: "chi2V",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_3_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_3_v,
            python: "chi_3_v", javascript: "chi3V",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_4_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_4_v,
            python: "chi_4_v", javascript: "chi4V",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_n_v",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_n_v,
            python: "chi_n_v", javascript: "chiNV",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: order, type: u32, default: required }],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, u32) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_0_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_0_n,
            python: "chi_0_n", javascript: "chi0N",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_1_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_1_n,
            python: "chi_1_n", javascript: "chi1N",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_2_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_2_n,
            python: "chi_2_n", javascript: "chi2N",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_3_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_3_n,
            python: "chi_3_n", javascript: "chi3N",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_4_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_4_n,
            python: "chi_4_n", javascript: "chi4N",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.chi_n_n",
            item: callable, owner: molecule,
            rust: crate::Molecule::chi_n_n,
            python: "chi_n_n", javascript: "chiNN",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: order, type: u32, default: required }],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, u32) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.molecular_weight",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::molecular_weight,
            python: "molecular_weight",
            javascript: "molecularWeight",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: f64,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.exact_molecular_weight",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::exact_molecular_weight,
            python: "exact_molecular_weight",
            javascript: "exactMolecularWeight",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: f64,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.molecular_formula",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::molecular_formula,
            python: "molecular_formula",
            javascript: "molecularFormula",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: String,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<String, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "types.RotatableBondsOptions", item: type, owner: type_,
            rust: crate::RotatableBondsOptions,
            python: "RotatableBondsOptions", javascript: "RotatableBondsOptions",
            feature: "cap-descriptors", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_amide_bonds", item: callable, owner: molecule,
            rust: crate::Molecule::num_amide_bonds, python: "num_amide_bonds", javascript: "numAmideBonds",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_spiro_atoms", item: callable, owner: molecule,
            rust: crate::Molecule::num_spiro_atoms, python: "num_spiro_atoms", javascript: "numSpiroAtoms",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_bridgehead_atoms", item: callable, owner: molecule,
            rust: crate::Molecule::num_bridgehead_atoms, python: "num_bridgehead_atoms", javascript: "numBridgeheadAtoms",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_atom_stereo_centers", item: callable, owner: molecule,
            rust: crate::Molecule::num_atom_stereo_centers, python: "num_atom_stereo_centers", javascript: "numAtomStereoCenters",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_unspecified_atom_stereo_centers", item: callable, owner: molecule,
            rust: crate::Molecule::num_unspecified_atom_stereo_centers, python: "num_unspecified_atom_stereo_centers", javascript: "numUnspecifiedAtomStereoCenters",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_rotatable_bonds", item: callable, owner: molecule,
            rust: crate::Molecule::num_rotatable_bonds, python: "num_rotatable_bonds", javascript: "numRotatableBonds",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_rotatable_bonds_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::num_rotatable_bonds_with_params, python: "num_rotatable_bonds_with_params", javascript: "numRotatableBondsWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: params, type: &crate::RotatableBondsOptions, default: required }],
            output: u32, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, &crate::RotatableBondsOptions) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.molecular_weight_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::molecular_weight_with_params, python: "molecular_weight_with_params", javascript: "molecularWeightWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: only_heavy, type: bool, default: required }],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.exact_molecular_weight_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::exact_molecular_weight_with_params, python: "exact_molecular_weight_with_params", javascript: "exactMolecularWeightWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: only_heavy, type: bool, default: required }],
            output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.molecular_formula_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::molecular_formula_with_params, python: "molecular_formula_with_params", javascript: "molecularFormulaWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: separate_isotopes, type: bool, default: required }, { name: abbreviate_h_isotopes, type: bool, default: required }],
            output: String, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool, bool) -> Result<String, crate::DescriptorReadError>,
        },
        // Descriptor count queries: the two public error types precede the
        // five registered query methods. DescriptorError is the canonical
        // re-export of the domain error; DescriptorReadError is the root
        // read boundary (missing prepared valence / typed algorithm cause).
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "types.DescriptorError",
            item: type,
            owner: type_,
            rust: crate::DescriptorError,
            python: "DescriptorError",
            javascript: "DescriptorError",
            feature: "cap-descriptors", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "types.DescriptorReadError",
            item: type,
            owner: type_,
            rust: crate::DescriptorReadError,
            python: "DescriptorReadError",
            javascript: "DescriptorReadError",
            feature: "cap-descriptors", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_heavy_atoms",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_heavy_atoms,
            python: "num_heavy_atoms",
            javascript: "numHeavyAtoms",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.total_atom_count",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::total_atom_count,
            python: "total_atom_count",
            javascript: "totalAtomCount",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_rings",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_rings,
            python: "num_rings",
            javascript: "numRings",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_heterocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_heterocycles,
            python: "num_heterocycles",
            javascript: "numHeterocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_heteroatoms",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_heteroatoms,
            python: "num_heteroatoms",
            javascript: "numHeteroatoms",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_hba",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_hba,
            python: "num_hba",
            javascript: "numHba",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_hbd",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_hbd,
            python: "num_hbd",
            javascript: "numHbd",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aromatic_rings",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aromatic_rings,
            python: "num_aromatic_rings",
            javascript: "numAromaticRings",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_saturated_rings",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_saturated_rings,
            python: "num_saturated_rings",
            javascript: "numSaturatedRings",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aliphatic_rings",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aliphatic_rings,
            python: "num_aliphatic_rings",
            javascript: "numAliphaticRings",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aromatic_heterocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aromatic_heterocycles,
            python: "num_aromatic_heterocycles",
            javascript: "numAromaticHeterocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aromatic_carbocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aromatic_carbocycles,
            python: "num_aromatic_carbocycles",
            javascript: "numAromaticCarbocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aliphatic_heterocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aliphatic_heterocycles,
            python: "num_aliphatic_heterocycles",
            javascript: "numAliphaticHeterocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_aliphatic_carbocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_aliphatic_carbocycles,
            python: "num_aliphatic_carbocycles",
            javascript: "numAliphaticCarbocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_saturated_heterocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_saturated_heterocycles,
            python: "num_saturated_heterocycles",
            javascript: "numSaturatedHeterocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.num_saturated_carbocycles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::num_saturated_carbocycles,
            python: "num_saturated_carbocycles",
            javascript: "numSaturatedCarbocycles",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.lipinski_hba",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::lipinski_hba,
            python: "lipinski_hba",
            javascript: "lipinskiHba",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.lipinski_hbd",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::lipinski_hbd,
            python: "lipinski_hbd",
            javascript: "lipinskiHbd",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: u32,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<u32, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.fraction_csp3",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::fraction_csp3,
            python: "fraction_csp3",
            javascript: "fractionCsp3",
            feature: "cap-descriptors", status: experimental,
            kind: instance,
            parameters: [],
            output: f64,
            error: crate::DescriptorReadError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.MmffProperties", item: type, owner: type_,
            rust: crate::MmffProperties, python: "MmffProperties", javascript: "MmffProperties",
            feature: "cap-forcefields", status: experimental, role: result,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.MmffPropertiesParams", item: type, owner: type_,
            rust: crate::MmffPropertiesParams, python: "MmffPropertiesParams", javascript: "MmffPropertiesParams",
            feature: "cap-forcefields", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.MmffAtomProperties", item: type, owner: type_,
            rust: crate::MmffAtomProperties, python: "MmffAtomProperties", javascript: "MmffAtomProperties",
            feature: "cap-forcefields", status: experimental, role: result,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.MmffVariant", item: type, owner: type_,
            rust: crate::MmffVariant, python: "MmffVariant", javascript: "MmffVariant",
            feature: "cap-forcefields", status: experimental, role: value,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.MmffMolPropertiesError", item: type, owner: type_,
            rust: crate::MmffMolPropertiesError, python: "MmffMolPropertiesError", javascript: "MmffMolPropertiesError",
            feature: "cap-forcefields", status: experimental, role: error,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.mmff_has_all_molecule_params", item: callable, owner: molecule,
            rust: crate::Molecule::mmff_has_all_molecule_params, python: "mmff_has_all_molecule_params", javascript: "mmffHasAllMoleculeParams",
            feature: "cap-forcefields", status: experimental, kind: instance,
            parameters: [], output: bool, error: crate::MmffMolPropertiesError,
            state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<bool, crate::MmffMolPropertiesError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.mmff_properties", item: callable, owner: molecule,
            rust: crate::Molecule::mmff_properties, python: "mmff_properties", javascript: "mmffProperties",
            feature: "cap-forcefields", status: experimental, kind: instance,
            parameters: [], output: crate::MmffProperties, error: crate::MmffMolPropertiesError,
            state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<crate::MmffProperties, crate::MmffMolPropertiesError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.mmff_properties_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::mmff_properties_with_params, python: "mmff_properties_with_params", javascript: "mmffPropertiesWithParams",
            feature: "cap-forcefields", status: experimental, kind: instance,
            parameters: [{ name: params, type: &crate::MmffPropertiesParams, default: required }], output: crate::MmffProperties, error: crate::MmffMolPropertiesError,
            state: read_only, operation: none, signature: fn(&crate::Molecule, &crate::MmffPropertiesParams) -> Result<crate::MmffProperties, crate::MmffMolPropertiesError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.uff_has_all_molecule_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::uff_has_all_molecule_params,
            python: "uff_has_all_molecule_params",
            javascript: "uffHasAllMoleculeParams",
            feature: "cap-forcefields",
            status: experimental,
            kind: instance,
            parameters: [],
            output: bool,
            error: crate::UffParameterQueryError,
            state: read_only,
            operation: none,
            signature: fn(&crate::Molecule) -> Result<bool, crate::UffParameterQueryError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffOptimizationParams", item: type, owner: type_,
            rust: crate::UffOptimizationParams, python: "UffOptimizationParams", javascript: "UffOptimizationParams",
            feature: "cap-forcefields", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffOptimizationResult", item: type, owner: type_,
            rust: crate::UffOptimizationResult, python: "UffOptimizationResult", javascript: "UffOptimizationResult",
            feature: "cap-forcefields", status: experimental, role: result,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffOptimizationError", item: type, owner: type_,
            rust: crate::UffOptimizationError, python: "UffOptimizationError", javascript: "UffOptimizationError",
            feature: "cap-forcefields", status: experimental, role: error,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffOptimizationErrorKind", item: type, owner: type_,
            rust: crate::UffOptimizationErrorKind, python: "UffOptimizationErrorKind", javascript: "UffOptimizationErrorKind",
            feature: "cap-forcefields", status: experimental, role: value,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffConformerOptimizationParams", item: type, owner: type_,
            rust: crate::UffConformerOptimizationParams, python: "UffConformerOptimizationParams", javascript: "UffConformerOptimizationParams",
            feature: "cap-forcefields", status: experimental, role: parameter,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffConformerOptimizationResult", item: type, owner: type_,
            rust: crate::UffConformerOptimizationResult, python: "UffConformerOptimizationResult", javascript: "UffConformerOptimizationResult",
            feature: "cap-forcefields", status: experimental, role: result,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "types.UffConformerResult", item: type, owner: type_,
            rust: crate::UffConformerResult, python: "UffConformerResult", javascript: "UffConformerResult",
            feature: "cap-forcefields", status: experimental, role: value,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "UffOptimizationError.kind", item: callable, owner: type_,
            rust: crate::UffOptimizationError::kind, python: "kind", javascript: "kind",
            feature: "cap-forcefields", status: experimental, kind: instance,
            parameters: [], output: crate::UffOptimizationErrorKind, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::UffOptimizationError) -> crate::UffOptimizationErrorKind,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.with_uff_optimized", item: callable, owner: molecule,
            rust: crate::Molecule::with_uff_optimized,
            python: "with_uff_optimized", javascript: "withUffOptimized",
            feature: "cap-forcefields", kind: instance,
            parameters: [], output: crate::UffOptimizationResult, error: crate::OperationError,
            state: value_returning, operation: "with_uff_optimized",
            signature: fn(&crate::Molecule) -> Result<crate::UffOptimizationResult, crate::OperationError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.with_uff_optimized_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::with_uff_optimized_with_params,
            python: "with_uff_optimized_with_params", javascript: "withUffOptimizedWithParams",
            feature: "cap-forcefields", kind: instance,
            parameters: [{ name: params, type: &crate::UffOptimizationParams, default: required }],
            output: crate::UffOptimizationResult, error: crate::OperationError,
            state: value_returning, operation: "with_uff_optimized_with_params",
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::UffOptimizationParams) -> Result<crate::UffOptimizationResult, crate::OperationError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.with_uff_optimized_confs", item: callable, owner: molecule,
            rust: crate::Molecule::with_uff_optimized_confs,
            python: "with_uff_optimized_confs", javascript: "withUffOptimizedConfs",
            feature: "cap-forcefields", kind: instance,
            parameters: [], output: crate::UffConformerOptimizationResult, error: crate::OperationError,
            state: value_returning, operation: "with_uff_optimized_confs",
            signature: fn(&crate::Molecule) -> Result<crate::UffConformerOptimizationResult, crate::OperationError>,
        },
        #[cfg(feature = "cap-forcefields")]
        {
            semantic_id: "Molecule.with_uff_optimized_confs_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::with_uff_optimized_confs_with_params,
            python: "with_uff_optimized_confs_with_params", javascript: "withUffOptimizedConfsWithParams",
            feature: "cap-forcefields", kind: instance,
            parameters: [{ name: params, type: &crate::UffConformerOptimizationParams, default: required }],
            output: crate::UffConformerOptimizationResult, error: crate::OperationError,
            state: value_returning, operation: "with_uff_optimized_confs_with_params",
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b crate::UffConformerOptimizationParams) -> Result<crate::UffConformerOptimizationResult, crate::OperationError>,
        },
        // H-add_stereo canonical public surface. The detached
        // AddHydrogensResult/HydrogenWarning hand-off remains internal; all
        // four callables project the same generated strong operation.
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "types.AddHsParams",
            item: type,
            owner: type_,
            rust: crate::AddHsParams,
            python: "AddHsParams",
            javascript: "AddHsParams",
            feature: "cap-hydrogens", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "types.HydrogenError",
            item: type,
            owner: type_,
            rust: crate::HydrogenError,
            python: "HydrogenError",
            javascript: "HydrogenError",
            feature: "cap-hydrogens", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "types.RemoveHsParams",
            item: type,
            owner: type_,
            rust: crate::RemoveHsParams,
            python: "RemoveHsParams",
            javascript: "RemoveHsParams",
            feature: "cap-hydrogens", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.with_hydrogens",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_hydrogens,
            python: "with_hydrogens",
            javascript: "withHydrogens",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_hydrogens",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.with_hydrogens_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::with_hydrogens_with_params,
            python: "with_hydrogens_with_params",
            javascript: "withHydrogensWithParams",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [
                { name: params, type: &crate::AddHsParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "with_hydrogens_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::AddHsParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.add_hydrogens_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::add_hydrogens_,
            python: "add_hydrogens_",
            javascript: "addHydrogens",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "add_hydrogens_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.add_hydrogens_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::add_hydrogens_with_params_,
            python: "add_hydrogens_with_params_",
            javascript: "addHydrogensWithParams",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [
                { name: params, type: &crate::AddHsParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "add_hydrogens_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::AddHsParams,
            ) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.without_hydrogens",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::without_hydrogens,
            python: "without_hydrogens",
            javascript: "withoutHydrogens",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "without_hydrogens",
            signature: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.without_hydrogens_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::without_hydrogens_with_params,
            python: "without_hydrogens_with_params",
            javascript: "withoutHydrogensWithParams",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [
                { name: params, type: &crate::RemoveHsParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::OperationError,
            state: value_returning,
            operation: "without_hydrogens_with_params",
            signature: for<'a, 'b> fn(
                &'a crate::Molecule,
                &'b crate::RemoveHsParams,
            ) -> Result<crate::Molecule, crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.remove_hydrogens_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::remove_hydrogens_,
            python: "remove_hydrogens_",
            javascript: "removeHydrogens",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "remove_hydrogens_",
            signature: fn(&mut crate::Molecule) -> Result<(), crate::OperationError>,
        },
        #[cfg(feature = "cap-hydrogens")]
        {
            semantic_id: "Molecule.remove_hydrogens_with_params_",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::remove_hydrogens_with_params_,
            python: "remove_hydrogens_with_params_",
            javascript: "removeHydrogensWithParams",
            feature: "cap-hydrogens",
            kind: instance,
            parameters: [
                { name: params, type: &crate::RemoveHsParams, default: required },
            ],
            output: (),
            error: crate::OperationError,
            state: in_place,
            operation: "remove_hydrogens_with_params_",
            signature: for<'a, 'b> fn(
                &'a mut crate::Molecule,
                &'b crate::RemoveHsParams,
            ) -> Result<(), crate::OperationError>,
        },
        // BIO-residue owns these immutable vocabulary/source-metadata values.
        // The facade only re-exports their unique detached owner.
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueInfoKind",
            item: type,
            owner: type_,
            rust: crate::ResidueInfoKind,
            python: "ResidueInfoKind",
            javascript: "ResidueInfoKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueInfoKind.name", item: callable,
            owner: type_, rust: crate::ResidueInfoKind::name,
            python: "name", javascript: "name",
            feature: "cap-bio", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: &'static str, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::ResidueInfoKind) -> &'static str,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueCode",
            item: type,
            owner: type_,
            rust: crate::ResidueCode,
            python: "ResidueCode",
            javascript: "ResidueCode",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueCodeParseError",
            item: type,
            owner: type_,
            rust: crate::ResidueCodeParseError,
            python: "ResidueCodeParseError",
            javascript: "ResidueCodeParseError",
            feature: "cap-bio", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueCodeParseError.input", item: callable,
            owner: type_, rust: crate::ResidueCodeParseError::input,
            python: "input", javascript: "input",
            feature: "cap-bio", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: &str, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::ResidueCodeParseError) -> &str,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueIdentity",
            item: type,
            owner: type_,
            rust: crate::ResidueIdentity,
            python: "ResidueIdentity",
            javascript: "ResidueIdentity",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueIdentity.new", item: callable,
            owner: type_, rust: crate::ResidueIdentity::new,
            python: "new", javascript: "new",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: name, type: String, default: required }],
            output: crate::ResidueIdentity, error: none,
            state: value_returning, operation: none,
            signature: fn(String) -> crate::ResidueIdentity,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueIdentity.name", item: callable,
            owner: type_, rust: crate::ResidueIdentity::name,
            python: "name", javascript: "name",
            feature: "cap-bio", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: &str, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::ResidueIdentity) -> &str,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueIdentity.code", item: callable,
            owner: type_, rust: crate::ResidueIdentity::code,
            python: "code", javascript: "code",
            feature: "cap-bio", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: crate::ResidueCode, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::ResidueIdentity) -> crate::ResidueCode,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueIdentity.info", item: callable,
            owner: type_, rust: crate::ResidueIdentity::info,
            python: "info", javascript: "info",
            feature: "cap-bio", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: crate::ResidueInfo, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::ResidueIdentity) -> crate::ResidueInfo,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueIdentity.is_tabulated", item: callable,
            owner: type_, rust: crate::ResidueIdentity::is_tabulated,
            python: "is_tabulated", javascript: "isTabulated",
            feature: "cap-bio", status: experimental, kind: instance, receiver: shared,
            parameters: [], output: bool, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::ResidueIdentity) -> bool,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueInfo",
            item: type,
            owner: type_,
            rust: crate::ResidueInfo,
            python: "ResidueInfo",
            javascript: "ResidueInfo",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueInfo.canonical_one_letter_code", item: callable,
            owner: type_, rust: crate::ResidueInfo::canonical_one_letter_code,
            python: "canonical_one_letter_code", javascript: "canonicalOneLetterCode",
            feature: "cap-bio", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: Option<char>, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::ResidueInfo) -> Option<char>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueInfo.parent_standard_code", item: callable,
            owner: type_, rust: crate::ResidueInfo::parent_standard_code,
            python: "parent_standard_code", javascript: "parentStandardCode",
            feature: "cap-bio", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: Option<crate::ResidueCode>, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::ResidueInfo) -> Option<crate::ResidueCode>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ResidueInfo.is_modified_amino_acid", item: callable,
            owner: type_, rust: crate::ResidueInfo::is_modified_amino_acid,
            python: "is_modified_amino_acid", javascript: "isModifiedAminoAcid",
            feature: "cap-bio", status: experimental, kind: instance, receiver: owned,
            parameters: [], output: bool, error: none,
            state: value_returning, operation: none,
            signature: fn(crate::ResidueInfo) -> bool,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.PdbAtomSerial",
            item: type,
            owner: type_,
            rust: crate::PdbAtomSerial,
            python: "PdbAtomSerial",
            javascript: "PdbAtomSerial",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.PdbChainId",
            item: type,
            owner: type_,
            rust: crate::PdbChainId,
            python: "PdbChainId",
            javascript: "PdbChainId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.PdbSeqId",
            item: type,
            owner: type_,
            rust: crate::PdbSeqId,
            python: "PdbSeqId",
            javascript: "PdbSeqId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.AtomName",
            item: type,
            owner: type_,
            rust: crate::AtomName,
            python: "AtomName",
            javascript: "AtomName",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueName",
            item: type,
            owner: type_,
            rust: crate::ResidueName,
            python: "ResidueName",
            javascript: "ResidueName",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.AltLocLabel",
            item: type,
            owner: type_,
            rust: crate::AltLocLabel,
            python: "AltLocLabel",
            javascript: "AltLocLabel",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.AtomSourceIds",
            item: type,
            owner: type_,
            rust: crate::AtomSourceIds,
            python: "AtomSourceIds",
            javascript: "AtomSourceIds",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueSourceIds",
            item: type,
            owner: type_,
            rust: crate::ResidueSourceIds,
            python: "ResidueSourceIds",
            javascript: "ResidueSourceIds",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ChainSourceIds",
            item: type,
            owner: type_,
            rust: crate::ChainSourceIds,
            python: "ChainSourceIds",
            javascript: "ChainSourceIds",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.EntitySourceIds",
            item: type,
            owner: type_,
            rust: crate::EntitySourceIds,
            python: "EntitySourceIds",
            javascript: "EntitySourceIds",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueSequenceError",
            item: type,
            owner: type_,
            rust: crate::ResidueSequenceError,
            python: "ResidueSequenceError",
            javascript: "ResidueSequenceError",
            feature: "cap-bio", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.residue_info",
            item: callable,
            owner: module,
            rust: crate::residue_info,
            python: "residue_info",
            javascript: "residueInfo",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [{ name: index, type: usize, default: required }],
            output: crate::ResidueInfo,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(usize) -> crate::ResidueInfo,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.residue_info_checked",
            item: callable,
            owner: module,
            rust: crate::residue_info_checked,
            python: "residue_info_checked",
            javascript: "residueInfoChecked",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [{ name: index, type: usize, default: required }],
            output: Option<crate::ResidueInfo>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(usize) -> Option<crate::ResidueInfo>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.find_residue_info_index",
            item: callable,
            owner: module,
            rust: crate::find_residue_info_index,
            python: "find_residue_info_index",
            javascript: "findResidueInfoIndex",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [{ name: name, type: &str, default: required }],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.find_residue_info",
            item: callable,
            owner: module,
            rust: crate::find_residue_info,
            python: "find_residue_info",
            javascript: "findResidueInfo",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [{ name: name, type: &str, default: required }],
            output: crate::ResidueInfo,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> crate::ResidueInfo,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.residue_code",
            item: callable,
            owner: module,
            rust: crate::residue_code,
            python: "residue_code",
            javascript: "residueCode",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [{ name: name, type: &str, default: required }],
            output: crate::ResidueCode,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&str) -> crate::ResidueCode,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.expand_one_letter",
            item: callable,
            owner: module,
            rust: crate::expand_one_letter,
            python: "expand_one_letter",
            javascript: "expandOneLetter",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [
                { name: code, type: char, default: required },
                { name: kind, type: crate::ResidueInfoKind, default: required },
            ],
            output: Option<&'static str>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(char, crate::ResidueInfoKind) -> Option<&'static str>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "module.expand_one_letter_sequence",
            item: callable,
            owner: module,
            rust: crate::expand_one_letter_sequence,
            python: "expand_one_letter_sequence",
            javascript: "expandOneLetterSequence",
            feature: "cap-bio", status: experimental,
            kind: module,
            parameters: [
                { name: sequence, type: &str, default: required },
                { name: kind, type: crate::ResidueInfoKind, default: required },
            ],
            output: Vec<String>,
            error: crate::ResidueSequenceError,
            state: read_only,
            operation: none,
            signature: fn(&str, crate::ResidueInfoKind) -> Result<Vec<String>, crate::ResidueSequenceError>,
        },
        // BioStructure is the public owner of detached BIO blocks. Algorithms
        // and structural validation remain in cosmolkit-bio; generated BIO
        // operations provide value/in-place wrappers over the same body.
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioStructure",
            item: type,
            owner: type_,
            rust: crate::BioStructure,
            python: "BioStructure",
            javascript: "BioStructure",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioStructureParts",
            item: type,
            owner: type_,
            rust: crate::BioStructureParts,
            python: "BioStructureParts",
            javascript: "BioStructureParts",
            feature: "cap-bio", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioStructureError",
            item: type,
            owner: type_,
            rust: crate::BioStructureError,
            python: "BioStructureError",
            javascript: "BioStructureError",
            feature: "cap-bio", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCoordinateFormat",
            item: type,
            owner: type_,
            rust: crate::BioCoordinateFormat,
            python: "BioCoordinateFormat",
            javascript: "BioCoordinateFormat",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCalcFlag",
            item: type,
            owner: type_,
            rust: crate::BioCalcFlag,
            python: "BioCalcFlag",
            javascript: "BioCalcFlag",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.EntityKind",
            item: type,
            owner: type_,
            rust: crate::EntityKind,
            python: "EntityKind",
            javascript: "EntityKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.PolymerKind",
            item: type,
            owner: type_,
            rust: crate::PolymerKind,
            python: "PolymerKind",
            javascript: "PolymerKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ResidueKind",
            item: type,
            owner: type_,
            rust: crate::ResidueKind,
            python: "ResidueKind",
            javascript: "ResidueKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ChainKind",
            item: type,
            owner: type_,
            rust: crate::ChainKind,
            python: "ChainKind",
            javascript: "ChainKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioRowSpan",
            item: type,
            owner: type_,
            rust: crate::BioRowSpan<crate::BioAtomId>,
            python: "BioRowSpan",
            javascript: "BioRowSpan",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAtomId",
            item: type,
            owner: type_,
            rust: crate::BioAtomId,
            python: "BioAtomId",
            javascript: "BioAtomId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioResidueId",
            item: type,
            owner: type_,
            rust: crate::BioResidueId,
            python: "BioResidueId",
            javascript: "BioResidueId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioChainId",
            item: type,
            owner: type_,
            rust: crate::BioChainId,
            python: "BioChainId",
            javascript: "BioChainId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioEntityId",
            item: type,
            owner: type_,
            rust: crate::BioEntityId,
            python: "BioEntityId",
            javascript: "BioEntityId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioModelId",
            item: type,
            owner: type_,
            rust: crate::BioModelId,
            python: "BioModelId",
            javascript: "BioModelId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAssemblyId",
            item: type,
            owner: type_,
            rust: crate::BioAssemblyId,
            python: "BioAssemblyId",
            javascript: "BioAssemblyId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAltLocGroupId",
            item: type,
            owner: type_,
            rust: crate::BioAltLocGroupId,
            python: "BioAltLocGroupId",
            javascript: "BioAltLocGroupId",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAtomRow",
            item: type,
            owner: type_,
            rust: crate::BioAtomRow,
            python: "BioAtomRow",
            javascript: "BioAtomRow",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioResidueRow",
            item: type,
            owner: type_,
            rust: crate::BioResidueRow,
            python: "BioResidueRow",
            javascript: "BioResidueRow",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioChainRow",
            item: type,
            owner: type_,
            rust: crate::BioChainRow,
            python: "BioChainRow",
            javascript: "BioChainRow",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioEntityRow",
            item: type,
            owner: type_,
            rust: crate::BioEntityRow,
            python: "BioEntityRow",
            javascript: "BioEntityRow",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioModelRow",
            item: type,
            owner: type_,
            rust: crate::BioModelRow,
            python: "BioModelRow",
            javascript: "BioModelRow",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCoordinateBlock",
            item: type,
            owner: type_,
            rust: crate::BioCoordinateBlock,
            python: "BioCoordinateBlock",
            javascript: "BioCoordinateBlock",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioTransform",
            item: type,
            owner: type_,
            rust: crate::BioTransform,
            python: "BioTransform",
            javascript: "BioTransform",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCrystalCell",
            item: type,
            owner: type_,
            rust: crate::BioCrystalCell,
            python: "BioCrystalCell",
            javascript: "BioCrystalCell",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCrystalInfo",
            item: type,
            owner: type_,
            rust: crate::BioCrystalInfo,
            python: "BioCrystalInfo",
            javascript: "BioCrystalInfo",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioNcsOperator",
            item: type,
            owner: type_,
            rust: crate::BioNcsOperator,
            python: "BioNcsOperator",
            javascript: "BioNcsOperator",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAssemblyOperator",
            item: type,
            owner: type_,
            rust: crate::BioAssemblyOperator,
            python: "BioAssemblyOperator",
            javascript: "BioAssemblyOperator",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAssemblyGenerator",
            item: type,
            owner: type_,
            rust: crate::BioAssemblyGenerator,
            python: "BioAssemblyGenerator",
            javascript: "BioAssemblyGenerator",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAssemblySpecialKind",
            item: type,
            owner: type_,
            rust: crate::BioAssemblySpecialKind,
            python: "BioAssemblySpecialKind",
            javascript: "BioAssemblySpecialKind",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioAssembly",
            item: type,
            owner: type_,
            rust: crate::BioAssembly,
            python: "BioAssembly",
            javascript: "BioAssembly",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.AltLocRequest",
            item: type,
            owner: type_,
            rust: crate::AltLocRequest,
            python: "AltLocRequest",
            javascript: "AltLocRequest",
            feature: "cap-bio", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioEntityDbRef",
            item: type,
            owner: type_,
            rust: crate::BioEntityDbRef,
            python: "BioEntityDbRef",
            javascript: "BioEntityDbRef",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSiftsUnpResidue",
            item: type,
            owner: type_,
            rust: crate::BioSiftsUnpResidue,
            python: "BioSiftsUnpResidue",
            javascript: "BioSiftsUnpResidue",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        // Protein uses the same BIO operation boundary while preserving its
        // amino-acid invariant. Child views borrow the detached hierarchy;
        // associated text constructors delegate parsing to cosmolkit-io.
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.Protein",
            item: type,
            owner: type_,
            rust: crate::Protein,
            python: "Protein",
            javascript: "Protein",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinProjectionError",
            item: type,
            owner: type_,
            rust: crate::ProteinProjectionError,
            python: "ProteinProjectionError",
            javascript: "ProteinProjectionError",
            feature: "cap-bio", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinSelectionSummary",
            item: type,
            owner: type_,
            rust: crate::ProteinSelectionSummary,
            python: "ProteinSelectionSummary",
            javascript: "ProteinSelectionSummary",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinChainIter", item: type, owner: type_,
            rust: crate::ProteinChainIter<'static>, python: "ProteinChainIter", javascript: "ProteinChainIter",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinResidueIter", item: type, owner: type_,
            rust: crate::ProteinResidueIter<'static>, python: "ProteinResidueIter", javascript: "ProteinResidueIter",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinAtomIter", item: type, owner: type_,
            rust: crate::ProteinAtomIter<'static>, python: "ProteinAtomIter", javascript: "ProteinAtomIter",
            feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinChainRef",
            item: type,
            owner: type_,
            rust: crate::ProteinChainRef<'static>,
            python: "ProteinChainRef",
            javascript: "ProteinChainRef",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinResidueRef",
            item: type,
            owner: type_,
            rust: crate::ProteinResidueRef<'static>,
            python: "ProteinResidueRef",
            javascript: "ProteinResidueRef",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.ProteinAtomRef",
            item: type,
            owner: type_,
            rust: crate::ProteinAtomRef<'static>,
            python: "ProteinAtomRef",
            javascript: "ProteinAtomRef",
            feature: "cap-bio", status: experimental,
            role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.from_parts", item: callable, owner: type_,
            rust: crate::BioStructure::from_parts, python: "from_parts", javascript: "fromParts",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: parts, type: crate::BioStructureParts, default: required }],
            output: crate::BioStructure, error: crate::BioStructureError, state: value_returning, operation: none,
            signature: fn(crate::BioStructureParts) -> Result<crate::BioStructure, crate::BioStructureError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.validate_parts", item: callable, owner: type_,
            rust: crate::BioStructure::validate_parts, python: "validate_parts", javascript: "validateParts",
            feature: "cap-bio", status: experimental, kind: static_,
            parameters: [{ name: parts, type: &crate::BioStructureParts, default: required }],
            output: (), error: crate::BioStructureError, state: read_only, operation: none,
            signature: for<'b> fn(&'b crate::BioStructureParts) -> Result<(), crate::BioStructureError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.validate", item: callable, owner: type_,
            rust: crate::BioStructure::validate, python: "validate", javascript: "validate",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: (), error: crate::BioStructureError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> Result<(), crate::BioStructureError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.into_parts", item: callable, owner: type_,
            rust: crate::BioStructure::into_parts, python: "into_parts", javascript: "intoParts",
            feature: "cap-bio", status: experimental, kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioStructureParts, error: none, state: value_returning, operation: none,
            signature: fn(crate::BioStructure) -> crate::BioStructureParts,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.input_format", item: callable, owner: type_,
            rust: crate::BioStructure::input_format, python: "input_format", javascript: "inputFormat",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: crate::BioCoordinateFormat, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> crate::BioCoordinateFormat,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.models", item: callable, owner: type_,
            rust: crate::BioStructure::models, python: "models", javascript: "models",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioModelRow], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioModelRow],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.chains", item: callable, owner: type_,
            rust: crate::BioStructure::chains, python: "chains", javascript: "chains",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioChainRow], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioChainRow],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.residues", item: callable, owner: type_,
            rust: crate::BioStructure::residues, python: "residues", javascript: "residues",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioResidueRow], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioResidueRow],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.atoms", item: callable, owner: type_,
            rust: crate::BioStructure::atoms, python: "atoms", javascript: "atoms",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioAtomRow], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioAtomRow],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.entities", item: callable, owner: type_,
            rust: crate::BioStructure::entities, python: "entities", javascript: "entities",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioEntityRow], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioEntityRow],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.connections", item: callable, owner: type_,
            rust: crate::BioStructure::connections, python: "connections", javascript: "connections",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioConnection], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioConnection],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.cispeps", item: callable, owner: type_,
            rust: crate::BioStructure::cispeps, python: "cispeps", javascript: "cispeps",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioCisPep], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioCisPep],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.mod_residues", item: callable, owner: type_,
            rust: crate::BioStructure::mod_residues, python: "mod_residues", javascript: "modResidues",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioModRes], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioModRes],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.helices", item: callable, owner: type_,
            rust: crate::BioStructure::helices, python: "helices", javascript: "helices",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioHelix], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioHelix],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.sheets", item: callable, owner: type_,
            rust: crate::BioStructure::sheets, python: "sheets", javascript: "sheets",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioSheet], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioSheet],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.ncs_operators", item: callable, owner: type_,
            rust: crate::BioStructure::ncs_operators, python: "ncs_operators", javascript: "ncsOperators",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioNcsOperator], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioNcsOperator],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.assemblies", item: callable, owner: type_,
            rust: crate::BioStructure::assemblies, python: "assemblies", javascript: "assemblies",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &[crate::BioAssembly], error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a [crate::BioAssembly],
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.metadata", item: callable, owner: type_,
            rust: crate::BioStructure::metadata, python: "metadata", javascript: "metadata",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &crate::BioMetadata, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a crate::BioMetadata,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.source_state", item: callable, owner: type_,
            rust: crate::BioStructure::source_state, python: "source_state", javascript: "sourceState",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &crate::BioStructureSourceState, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a crate::BioStructureSourceState,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.coordinates", item: callable, owner: type_,
            rust: crate::BioStructure::coordinates, python: "coordinates", javascript: "coordinates",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: &crate::BioCoordinateBlock, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a crate::BioCoordinateBlock,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.crystal", item: callable, owner: type_,
            rust: crate::BioStructure::crystal, python: "crystal", javascript: "crystal",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [],
            output: Option<&crate::BioCrystalInfo>, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> Option<&'a crate::BioCrystalInfo>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.find_entity", item: callable, owner: type_,
            rust: crate::BioStructure::find_entity, python: "find_entity", javascript: "findEntity",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: source_id, type: &str, default: required }],
            output: Option<(crate::BioEntityId, &crate::BioEntityRow)>, error: none, state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::BioStructure, &'b str) -> Option<(crate::BioEntityId, &'a crate::BioEntityRow)>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.find_entity_of_subchain", item: callable, owner: type_,
            rust: crate::BioStructure::find_entity_of_subchain, python: "find_entity_of_subchain", javascript: "findEntityOfSubchain",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: subchain, type: &str, default: required }],
            output: Option<(crate::BioEntityId, &crate::BioEntityRow)>, error: none, state: read_only, operation: none,
            signature: for<'a, 'b> fn(&'a crate::BioStructure, &'b str) -> Option<(crate::BioEntityId, &'a crate::BioEntityRow)>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.find_atom", item: callable, owner: type_,
            rust: crate::BioStructure::find_atom, python: "find_atom", javascript: "findAtom",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: residue_id, type: crate::BioResidueId, default: required }, { name: name, type: crate::AtomName, default: required }, { name: request, type: crate::AltLocRequest, default: required }, { name: element, type: Option<crate::Element>, default: required }],
            output: Option<(crate::BioAtomId, &crate::BioAtomRow)>, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, crate::BioResidueId, crate::AtomName, crate::AltLocRequest, Option<crate::Element>) -> Option<(crate::BioAtomId, &'a crate::BioAtomRow)>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.atom_by_altloc", item: callable, owner: type_,
            rust: crate::BioStructure::atom_by_altloc, python: "atom_by_altloc", javascript: "atomByAltloc",
            feature: "cap-bio", status: experimental, kind: instance,
            parameters: [{ name: residue_id, type: crate::BioResidueId, default: required }, { name: name, type: crate::AtomName, default: required }, { name: altloc, type: Option<crate::AltLocLabel>, default: required }],
            output: (crate::BioAtomId, &crate::BioAtomRow), error: crate::BioStructureError, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, crate::BioResidueId, crate::AtomName, Option<crate::AltLocLabel>) -> Result<(crate::BioAtomId, &'a crate::BioAtomRow), crate::BioStructureError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioConnection", item: type, owner: type_, rust: crate::BioConnection,
            python: "BioConnection", javascript: "BioConnection", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioCisPep", item: type, owner: type_, rust: crate::BioCisPep,
            python: "BioCisPep", javascript: "BioCisPep", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioModRes", item: type, owner: type_, rust: crate::BioModRes,
            python: "BioModRes", javascript: "BioModRes", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioHelix", item: type, owner: type_, rust: crate::BioHelix,
            python: "BioHelix", javascript: "BioHelix", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioSheet", item: type, owner: type_, rust: crate::BioSheet,
            python: "BioSheet", javascript: "BioSheet", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioMetadata", item: type, owner: type_, rust: crate::BioMetadata,
            python: "BioMetadata", javascript: "BioMetadata", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "types.BioStructureSourceState", item: type, owner: type_, rust: crate::BioStructureSourceState,
            python: "BioStructureSourceState", javascript: "BioStructureSourceState", feature: "cap-bio", status: experimental, role: value,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.num_models",
            item: callable, owner: type_, rust: crate::BioStructure::num_models,
            python: "num_models", javascript: "numModels", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: usize, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.name", item: callable, owner: type_, rust: crate::BioStructure::name,
            python: "name", javascript: "name", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: &str, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a str,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.has_origx", item: callable, owner: type_, rust: crate::BioStructure::has_origx,
            python: "has_origx", javascript: "hasOrigx", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: bool, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> bool,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.origx", item: callable, owner: type_, rust: crate::BioStructure::origx,
            python: "origx", javascript: "origx", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: &crate::BioTransform, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> &'a crate::BioTransform,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.ncs_oper_identity_id", item: callable, owner: type_, rust: crate::BioStructure::ncs_oper_identity_id,
            python: "ncs_oper_identity_id", javascript: "ncsOperIdentityId", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: Option<&str>, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure) -> Option<&'a str>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.resolution", item: callable, owner: type_, rust: crate::BioStructure::resolution,
            python: "resolution", javascript: "resolution", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: f64, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> f64,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.ter_status", item: callable, owner: type_, rust: crate::BioStructure::ter_status,
            python: "ter_status", javascript: "terStatus", feature: "cap-bio", status: experimental, kind: instance,
            parameters: [], output: u8, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> u8,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.num_chains",
            item: callable, owner: type_, rust: crate::BioStructure::num_chains,
            python: "num_chains", javascript: "numChains", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: usize, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.num_residues",
            item: callable, owner: type_, rust: crate::BioStructure::num_residues,
            python: "num_residues", javascript: "numResidues", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: usize, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.num_atoms",
            item: callable, owner: type_, rust: crate::BioStructure::num_atoms,
            python: "num_atoms", javascript: "numAtoms", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: usize, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.num_entities",
            item: callable, owner: type_, rust: crate::BioStructure::num_entities,
            python: "num_entities", javascript: "numEntities", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: usize, error: none,
            state: read_only, operation: none,
            signature: fn(&crate::BioStructure) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.atom_position",
            item: callable, owner: type_, rust: crate::BioStructure::atom_position,
            python: "atom_position", javascript: "atomPosition", feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [{ name: atom, type: crate::BioAtomId, default: required }],
            output: Option<[f64; 3]>, error: none, state: read_only, operation: none,
            signature: fn(&crate::BioStructure, crate::BioAtomId) -> Option<[f64; 3]>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.residue_atoms",
            item: callable, owner: type_, rust: crate::BioStructure::residue_atoms,
            python: "residue_atoms", javascript: "residueAtoms", feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [{ name: residue, type: crate::BioResidueId, default: required }],
            output: Option<&[crate::BioAtomRow]>, error: none, state: read_only, operation: none,
            signature: for<'a> fn(&'a crate::BioStructure, crate::BioResidueId) -> Option<&'a [crate::BioAtomRow]>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "BioStructure.protein",
            item: callable,
            owner: type_,
            rust: crate::BioStructure::protein,
            python: "protein",
            javascript: "protein",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::Protein,
            error: crate::ProteinProjectionError,
            state: read_only,
            operation: none,
            signature: fn(&crate::BioStructure) -> Result<crate::Protein, crate::ProteinProjectionError>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.input_format",
            item: callable, owner: type_, rust: crate::Protein::input_format,
            python: "input_format", javascript: "inputFormat", feature: "cap-bio", status: experimental,
            kind: instance, parameters: [], output: crate::BioCoordinateFormat,
            error: none, state: read_only, operation: none,
            signature: fn(&crate::Protein) -> crate::BioCoordinateFormat,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.as_bio_structure",
            item: callable,
            owner: type_,
            rust: crate::Protein::as_bio_structure,
            python: "as_bio_structure",
            javascript: "asBioStructure",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::BioStructure,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::Protein) -> &'a crate::BioStructure,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.into_bio_structure",
            item: callable,
            owner: type_,
            rust: crate::Protein::into_bio_structure,
            python: "into_bio_structure",
            javascript: "intoBioStructure",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioStructure,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::Protein) -> crate::BioStructure,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.selection_summary",
            item: callable,
            owner: type_,
            rust: crate::Protein::selection_summary,
            python: "selection_summary",
            javascript: "selectionSummary",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::ProteinSelectionSummary,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein) -> crate::ProteinSelectionSummary,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.num_models",
            item: callable,
            owner: type_,
            rust: crate::Protein::num_models,
            python: "num_models",
            javascript: "numModels",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.num_chains",
            item: callable,
            owner: type_,
            rust: crate::Protein::num_chains,
            python: "num_chains",
            javascript: "numChains",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.num_residues",
            item: callable,
            owner: type_,
            rust: crate::Protein::num_residues,
            python: "num_residues",
            javascript: "numResidues",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.num_atoms",
            item: callable,
            owner: type_,
            rust: crate::Protein::num_atoms,
            python: "num_atoms",
            javascript: "numAtoms",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: usize,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein) -> usize,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.chains",
            item: callable,
            owner: type_,
            rust: crate::Protein::chains,
            python: "chains",
            javascript: "chains",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::ProteinChainIter<'_>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::Protein) -> crate::ProteinChainIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.chain",
            item: callable,
            owner: type_,
            rust: crate::Protein::chain,
            python: "chain",
            javascript: "chain",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [
                { name: index, type: usize, default: required },
            ],
            output: Option<crate::ProteinChainRef<'_>>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::Protein, usize) -> Option<crate::ProteinChainRef<'_>>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.residues",
            item: callable,
            owner: type_,
            rust: crate::Protein::residues,
            python: "residues",
            javascript: "residues",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::ProteinResidueIter<'_>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::Protein) -> crate::ProteinResidueIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "Protein.atoms",
            item: callable,
            owner: type_,
            rust: crate::Protein::atoms,
            python: "atoms",
            javascript: "atoms",
            feature: "cap-bio", status: experimental,
            kind: instance,
            parameters: [],
            output: crate::ProteinAtomIter<'_>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a> fn(&'a crate::Protein) -> crate::ProteinAtomIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.id",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::id,
            python: "id",
            javascript: "id",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioChainId,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinChainRef<'static>) -> crate::BioChainId,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.row",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::row,
            python: "row",
            javascript: "row",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: &crate::BioChainRow,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinChainRef<'a>) -> &'a crate::BioChainRow,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.kind",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::kind,
            python: "kind",
            javascript: "kind",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ChainKind,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinChainRef<'static>) -> crate::ChainKind,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.source",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::source,
            python: "source",
            javascript: "source",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: &crate::ChainSourceIds,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinChainRef<'a>) -> &'a crate::ChainSourceIds,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.residues",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::residues,
            python: "residues",
            javascript: "residues",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ProteinResidueIter<'_>,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinChainRef<'a>) -> crate::ProteinResidueIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinChainRef.atoms",
            item: callable,
            owner: type_,
            rust: crate::ProteinChainRef::atoms,
            python: "atoms",
            javascript: "atoms",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ProteinAtomIter<'_>,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinChainRef<'a>) -> crate::ProteinAtomIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.id",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::id,
            python: "id",
            javascript: "id",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioResidueId,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> crate::BioResidueId,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.row",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::row,
            python: "row",
            javascript: "row",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: &crate::BioResidueRow,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinResidueRef<'a>) -> &'a crate::BioResidueRow,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.name",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::name,
            python: "name",
            javascript: "name",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ResidueName,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> crate::ResidueName,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.kind",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::kind,
            python: "kind",
            javascript: "kind",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ResidueKind,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> crate::ResidueKind,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.info",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::info,
            python: "info",
            javascript: "info",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ResidueInfo,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> crate::ResidueInfo,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.code",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::code,
            python: "code",
            javascript: "code",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ResidueCode,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> crate::ResidueCode,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.one_letter_code",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::one_letter_code,
            python: "one_letter_code",
            javascript: "oneLetterCode",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: char,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> char,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.fasta_code",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::fasta_code,
            python: "fasta_code",
            javascript: "fastaCode",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: char,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> char,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.is_standard",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::is_standard,
            python: "is_standard",
            javascript: "isStandard",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: bool,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinResidueRef<'static>) -> bool,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.chain",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::chain,
            python: "chain",
            javascript: "chain",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ProteinChainRef<'_>,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinResidueRef<'a>) -> crate::ProteinChainRef<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinResidueRef.atoms",
            item: callable,
            owner: type_,
            rust: crate::ProteinResidueRef::atoms,
            python: "atoms",
            javascript: "atoms",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ProteinAtomIter<'_>,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinResidueRef<'a>) -> crate::ProteinAtomIter<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.id",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::id,
            python: "id",
            javascript: "id",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioAtomId,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinAtomRef<'static>) -> crate::BioAtomId,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.row",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::row,
            python: "row",
            javascript: "row",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: &crate::BioAtomRow,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinAtomRef<'a>) -> &'a crate::BioAtomRow,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.name",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::name,
            python: "name",
            javascript: "name",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::AtomName,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinAtomRef<'static>) -> crate::AtomName,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.element",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::element,
            python: "element",
            javascript: "element",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::Element,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinAtomRef<'static>) -> crate::Element,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.altloc",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::altloc,
            python: "altloc",
            javascript: "altloc",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: Option<crate::AltLocLabel>,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinAtomRef<'static>) -> Option<crate::AltLocLabel>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.residue",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::residue,
            python: "residue",
            javascript: "residue",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::ProteinResidueRef<'_>,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::ProteinAtomRef<'a>) -> crate::ProteinResidueRef<'a>,
        },
        #[cfg(feature = "cap-bio")]
        {
            semantic_id: "ProteinAtomRef.position",
            item: callable,
            owner: type_,
            rust: crate::ProteinAtomRef::position,
            python: "position",
            javascript: "position",
            feature: "cap-bio", status: experimental,
            kind: instance, receiver: owned,
            parameters: [],
            output: [f64; 3],
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::ProteinAtomRef<'static>) -> [f64; 3],
        },
        // Canonical Rust SDF construction contracts. The public SdfReadParams
        // contract specifies sanitize=true,
        // remove_hydrogens=true, expand_attachment_points=false,
        // process_property_lists=true, and coordinate_mode=Preserve.
        //
        // Molecule readers are concrete-only: after the unique detached
        // finalization pipeline, query-bearing results produce a structured
        // SdfError query-record category, never a lossy Molecule conversion.
        // SdfRecord readers preserve either Molecule or QueryGraph in SdfGraph,
        // together with ordered data fields, properties, SGroups and source
        // coordinate dimension. SdfGraph is an IO payload tag, not a third
        // chemistry model. Classification is made after finalization, including
        // source-defined attachment promotion; a no-op option is not promotion.
        // Only runtime constructs live Molecule; domain IO remains detached.
        // The SDF readers below are implemented in the public Rust crate.
        // CK-COORD-001 is the approved terminal-coordinate exception to exact
        // RDKit parity; the remaining pinned source behavior is still required.
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.SdfRecord",
            item: type,
            owner: type_,
            rust: crate::SdfRecord,
            python: "SdfRecord",
            javascript: "SdfRecord",
            feature: "cap-io", status: experimental,
            role: result,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.SdfGraph",
            item: type,
            owner: type_,
            rust: crate::SdfGraph,
            python: "SdfGraph",
            javascript: "SdfGraph",
            feature: "cap-io", status: experimental,
            role: result,
        },
        // The record payload is classified only after molfile finalization.
        // `molecule`/`query_graph` expose checked borrows and project a
        // structured SdfError::WrongGraphKind (expected/actual payload kind).
        // Concrete-only Molecule readers instead project SdfError::QueryRecord.
        // Every registered row resolves to its actual public Rust item.
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.graph",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::graph,
            python: "graph",
            javascript: "graph",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::SdfGraph,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> &crate::SdfGraph,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.molecule",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::molecule,
            python: "molecule",
            javascript: "molecule",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::Molecule,
            error: crate::SdfError,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> Result<&crate::Molecule, crate::SdfError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.query_graph",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::query_graph,
            python: "query_graph",
            javascript: "queryGraph",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::QueryGraph,
            error: crate::SdfError,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> Result<&crate::QueryGraph, crate::SdfError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.data_fields",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::data_fields,
            python: "data_fields",
            javascript: "dataFields",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &[(String, String)],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> &[(String, String)],
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.properties",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::properties,
            python: "properties",
            javascript: "properties",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &crate::MoleculeProperties,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> &crate::MoleculeProperties,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.substance_groups",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::substance_groups,
            python: "substance_groups",
            javascript: "substanceGroups",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: &[crate::SubstanceGroup],
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> &[crate::SubstanceGroup],
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.source_coordinate_dim",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::source_coordinate_dim,
            python: "source_coordinate_dim",
            javascript: "sourceCoordinateDim",
            feature: "cap-io", status: experimental,
            kind: instance,
            parameters: [],
            output: Option<crate::CoordinateDimension>,
            error: none,
            state: read_only,
            operation: none,
            signature: fn(&crate::SdfRecord) -> Option<crate::CoordinateDimension>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.from_sdf",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::from_sdf,
            python: "from_sdf",
            javascript: "fromSdf",
            feature: "cap-io", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
            ],
            output: crate::SdfRecord,
            error: crate::SdfError,
            state: value_returning,
            operation: none,
            signature: fn(&str) -> Result<crate::SdfRecord, crate::SdfError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "SdfRecord.from_sdf_with_params",
            item: callable,
            owner: type_,
            rust: crate::SdfRecord::from_sdf_with_params,
            python: "from_sdf_with_params",
            javascript: "fromSdfWithParams",
            feature: "cap-io", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
                { name: params, type: &crate::SdfReadParams, default: required },
            ],
            output: crate::SdfRecord,
            error: crate::SdfError,
            state: value_returning,
            operation: none,
            signature: for<'a, 'b> fn(&'a str, &'b crate::SdfReadParams) -> Result<crate::SdfRecord, crate::SdfError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.SdfCoordinateMode",
            item: type,
            owner: type_,
            rust: crate::SdfCoordinateMode,
            python: "SdfCoordinateMode",
            javascript: "SdfCoordinateMode",
            feature: "cap-io", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.SdfReadParams",
            item: type,
            owner: type_,
            rust: crate::SdfReadParams,
            python: "SdfReadParams",
            javascript: "SdfReadParams",
            feature: "cap-io", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "types.SdfError",
            item: type,
            owner: type_,
            rust: crate::SdfError,
            python: "SdfError",
            javascript: "SdfError",
            feature: "cap-io", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "Molecule.from_sdf",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::from_sdf,
            python: "from_sdf",
            javascript: "fromSdf",
            feature: "cap-io", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
            ],
            output: crate::Molecule,
            error: crate::SdfError,
            state: value_returning,
            operation: none,
            signature: fn(&str) -> Result<crate::Molecule, crate::SdfError>,
        },
        #[cfg(feature = "cap-io")]
        {
            semantic_id: "Molecule.from_sdf_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::from_sdf_with_params,
            python: "from_sdf_with_params",
            javascript: "fromSdfWithParams",
            feature: "cap-io", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
                { name: params, type: &crate::SdfReadParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::SdfError,
            state: value_returning,
            operation: none,
            signature: for<'a, 'b> fn(
                &'a str,
                &'b crate::SdfReadParams,
            ) -> Result<crate::Molecule, crate::SdfError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.SmilesParseParams",
            item: type,
            owner: type_,
            rust: crate::SmilesParseParams,
            python: "SmilesParseParams",
            javascript: "SmilesParseParams",
            feature: "cap-smiles", status: experimental,
            role: parameter,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.SmilesError",
            item: type,
            owner: type_,
            rust: crate::SmilesError,
            python: "SmilesError",
            javascript: "SmilesError",
            feature: "cap-smiles", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "types.SmilesStereoError",
            item: type,
            owner: type_,
            rust: crate::SmilesStereoError,
            python: "SmilesStereoError",
            javascript: "SmilesStereoError",
            feature: "cap-smiles", status: experimental,
            role: error,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.from_smiles",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::from_smiles,
            python: "from_smiles",
            javascript: "fromSmiles",
            feature: "cap-smiles", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
            ],
            output: crate::Molecule,
            error: crate::SmilesError,
            state: value_returning,
            operation: none,
            signature: fn(&str) -> Result<crate::Molecule, crate::SmilesError>,
        },
        #[cfg(feature = "cap-smiles")]
        {
            semantic_id: "Molecule.from_smiles_with_params",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::from_smiles_with_params,
            python: "from_smiles_with_params",
            javascript: "fromSmilesWithParams",
            feature: "cap-smiles", status: experimental,
            kind: static_,
            parameters: [
                { name: input, type: &str, default: required },
                { name: params, type: &crate::SmilesParseParams, default: required },
            ],
            output: crate::Molecule,
            error: crate::SmilesError,
            state: value_returning,
            operation: none,
            signature: for<'a, 'b> fn(
                &'a str,
                &'b crate::SmilesParseParams,
            ) -> Result<crate::Molecule, crate::SmilesError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "types.CrippenTotals", item: type, owner: type_, rust: crate::CrippenTotals,
            python: "CrippenTotals", javascript: "CrippenTotals", feature: "cap-descriptors", status: experimental, role: result,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "types.LabuteAsaContributions", item: type, owner: type_, rust: crate::LabuteAsaContributions,
            python: "LabuteAsaContributions", javascript: "LabuteAsaContributions", feature: "cap-descriptors", status: experimental, role: result,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.crippen_descriptors", item: callable, owner: molecule,
            rust: crate::Molecule::crippen_descriptors, python: "crippen_descriptors", javascript: "crippenDescriptors",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: crate::CrippenTotals, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<crate::CrippenTotals, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.labute_asa", item: callable, owner: molecule,
            rust: crate::Molecule::labute_asa, python: "labute_asa", javascript: "labuteAsa",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.labute_asa_contributions", item: callable, owner: molecule,
            rust: crate::Molecule::labute_asa_contributions, python: "labute_asa_contributions", javascript: "labuteAsaContributions",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: crate::LabuteAsaContributions, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<crate::LabuteAsaContributions, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.tpsa", item: callable, owner: molecule,
            rust: crate::Molecule::tpsa, python: "tpsa", javascript: "tpsa",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa, python: "slogp_vsa", javascript: "slogpVsa",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: Vec<f64>, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<Vec<f64>, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa, python: "smr_vsa", javascript: "smrVsa",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: Vec<f64>, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<Vec<f64>, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_1", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_1, python: "slogp_vsa_1", javascript: "slogpVsa1",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_2", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_2, python: "slogp_vsa_2", javascript: "slogpVsa2",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_3", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_3, python: "slogp_vsa_3", javascript: "slogpVsa3",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_4", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_4, python: "slogp_vsa_4", javascript: "slogpVsa4",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_5", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_5, python: "slogp_vsa_5", javascript: "slogpVsa5",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_6", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_6, python: "slogp_vsa_6", javascript: "slogpVsa6",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_7", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_7, python: "slogp_vsa_7", javascript: "slogpVsa7",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_8", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_8, python: "slogp_vsa_8", javascript: "slogpVsa8",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_9", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_9, python: "slogp_vsa_9", javascript: "slogpVsa9",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_10", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_10, python: "slogp_vsa_10", javascript: "slogpVsa10",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_11", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_11, python: "slogp_vsa_11", javascript: "slogpVsa11",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_12", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_12, python: "slogp_vsa_12", javascript: "slogpVsa12",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_1", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_1, python: "smr_vsa_1", javascript: "smrVsa1",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_2", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_2, python: "smr_vsa_2", javascript: "smrVsa2",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_3", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_3, python: "smr_vsa_3", javascript: "smrVsa3",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_4", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_4, python: "smr_vsa_4", javascript: "smrVsa4",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_5", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_5, python: "smr_vsa_5", javascript: "smrVsa5",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_6", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_6, python: "smr_vsa_6", javascript: "smrVsa6",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_7", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_7, python: "smr_vsa_7", javascript: "smrVsa7",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_8", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_8, python: "smr_vsa_8", javascript: "smrVsa8",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_9", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_9, python: "smr_vsa_9", javascript: "smrVsa9",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_10", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_10, python: "smr_vsa_10", javascript: "smrVsa10",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.crippen_descriptors_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::crippen_descriptors_with_params, python: "crippen_descriptors_with_params", javascript: "crippenDescriptorsWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: include_hydrogens, type: bool, default: required }, { name: force, type: bool, default: required }], output: crate::CrippenTotals, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool, bool) -> Result<crate::CrippenTotals, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.labute_asa_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::labute_asa_with_params, python: "labute_asa_with_params", javascript: "labuteAsaWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: include_hydrogens, type: bool, default: required }, { name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool, bool) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.labute_asa_contributions_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::labute_asa_contributions_with_params, python: "labute_asa_contributions_with_params", javascript: "labuteAsaContributionsWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: include_hydrogens, type: bool, default: required }, { name: force, type: bool, default: required }], output: crate::LabuteAsaContributions, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool, bool) -> Result<crate::LabuteAsaContributions, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.tpsa_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::tpsa_with_params, python: "tpsa_with_params", javascript: "tpsaWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: include_sulfur_phosphorus, type: bool, default: required }, { name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, bool, bool) -> Result<f64, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.slogp_vsa_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::slogp_vsa_with_params, python: "slogp_vsa_with_params", javascript: "slogpVsaWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: bins, type: Option<&[f64]>, default: required }, { name: force, type: bool, default: required }], output: Vec<f64>, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, Option<&[f64]>, bool) -> Result<Vec<f64>, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        {
            semantic_id: "Molecule.smr_vsa_with_params", item: callable, owner: molecule,
            rust: crate::Molecule::smr_vsa_with_params, python: "smr_vsa_with_params", javascript: "smrVsaWithParams",
            feature: "cap-descriptors", status: experimental,
            kind: instance, parameters: [{ name: bins, type: Option<&[f64]>, default: required }, { name: force, type: bool, default: required }], output: Vec<f64>, error: crate::DescriptorReadError,
            state: read_only, operation: none,
            signature: fn(&crate::Molecule, Option<&[f64]>, bool) -> Result<Vec<f64>, crate::DescriptorReadError>,
        },
        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.qed", item: callable, owner: molecule,
          rust: crate::Molecule::qed, python: "qed", javascript: "qed",
          feature: "cap-descriptors", status: experimental, kind: instance, parameters: [],
          output: f64, error: crate::DescriptorReadError, state: read_only, operation: none,
          signature: fn(&crate::Molecule) -> Result<f64, crate::DescriptorReadError>, },
        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_0_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_0_v_with_params, python: "chi_0_v_with_params", javascript: "chi0VWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_1_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_1_v_with_params, python: "chi_1_v_with_params", javascript: "chi1VWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_2_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_2_v_with_params, python: "chi_2_v_with_params", javascript: "chi2VWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_3_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_3_v_with_params, python: "chi_3_v_with_params", javascript: "chi3VWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_4_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_4_v_with_params, python: "chi_4_v_with_params", javascript: "chi4VWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_n_v_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_n_v_with_params, python: "chi_n_v_with_params", javascript: "chiNVWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: order, type: u32, default: required }, { name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, u32, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_0_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_0_n_with_params, python: "chi_0_n_with_params", javascript: "chi0NWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_1_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_1_n_with_params, python: "chi_1_n_with_params", javascript: "chi1NWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_2_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_2_n_with_params, python: "chi_2_n_with_params", javascript: "chi2NWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_3_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_3_n_with_params, python: "chi_3_n_with_params", javascript: "chi3NWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_4_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_4_n_with_params, python: "chi_4_n_with_params", javascript: "chi4NWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, bool) -> Result<f64, crate::DescriptorReadError>, },

        #[cfg(feature = "cap-descriptors")]
        { semantic_id: "Molecule.chi_n_n_with_params", item: callable, owner: molecule,
          rust: crate::Molecule::chi_n_n_with_params, python: "chi_n_n_with_params", javascript: "chiNNWithParams",
          feature: "cap-descriptors", status: experimental, kind: instance,
          parameters: [{ name: order, type: u32, default: required }, { name: force, type: bool, default: required }], output: f64, error: crate::DescriptorReadError,
          state: read_only, operation: none,
          signature: fn(&crate::Molecule, u32, bool) -> Result<f64, crate::DescriptorReadError>, },
        { semantic_id:"types.Atom",item:type,owner:type_,rust:crate::Atom,python:"Atom",javascript:"Atom",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.Bond",item:type,owner:type_,rust:crate::Bond,python:"Bond",javascript:"Bond",feature:"metadata",status:experimental,role:value, },
        {semantic_id:"types.Hybridization",item:type,owner:type_,rust:crate::Hybridization,python:"Hybridization",javascript:"Hybridization",feature:"metadata",status:experimental,role:value,},
        {semantic_id:"Atom.hybridization",item:callable,owner:type_,rust:crate::Atom::hybridization,python:"hybridization",javascript:"hybridization",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::Hybridization,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->crate::Hybridization,},
        { semantic_id:"types.BondOrder",item:type,owner:type_,rust:crate::BondOrder,python:"BondOrder",javascript:"BondOrder",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.BondDirection",item:type,owner:type_,rust:crate::BondDirection,python:"BondDirection",javascript:"BondDirection",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.BondStereo",item:type,owner:type_,rust:crate::BondStereo,python:"BondStereo",javascript:"BondStereo",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.ChiralTag",item:type,owner:type_,rust:crate::ChiralTag,python:"ChiralTag",javascript:"ChiralTag",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.AtomSpec",item:type,owner:type_,rust:crate::AtomSpec,python:"AtomSpec",javascript:"AtomSpec",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"types.BondSpec",item:type,owner:type_,rust:crate::BondSpec,python:"BondSpec",javascript:"BondSpec",feature:"metadata",status:experimental,role:value, },
        { semantic_id:"Atom.id",item:callable,owner:type_,rust:crate::Atom::id,python:"id",javascript:"id",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::AtomId,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->crate::AtomId, },
        { semantic_id:"Atom.element",item:callable,owner:type_,rust:crate::Atom::element,python:"element",javascript:"element",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::Element,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->crate::Element, },
        { semantic_id:"Atom.atomic_number",item:callable,owner:type_,rust:crate::Atom::atomic_number,python:"atomic_number",javascript:"atomicNumber",feature:"runtime",status:experimental,kind:instance,parameters:[],output:u8,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->u8, },
        { semantic_id:"Atom.formal_charge",item:callable,owner:type_,rust:crate::Atom::formal_charge,python:"formal_charge",javascript:"formalCharge",feature:"runtime",status:experimental,kind:instance,parameters:[],output:i8,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->i8, },
        { semantic_id:"Atom.chiral_tag",item:callable,owner:type_,rust:crate::Atom::chiral_tag,python:"chiral_tag",javascript:"chiralTag",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::ChiralTag,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->crate::ChiralTag, },
        { semantic_id:"Atom.chiral_tag_code",item:callable,owner:type_,rust:crate::Atom::chiral_tag_code,python:"chiral_tag_code",javascript:"chiralTagCode",feature:"runtime",status:experimental,kind:instance,parameters:[],output:i64,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->i64, },
        { semantic_id:"Atom.chiral_tag_name",item:callable,owner:type_,rust:crate::Atom::chiral_tag_name,python:"chiral_tag_name",javascript:"chiralTagName",feature:"runtime",status:experimental,kind:instance,parameters:[],output:&'static str,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->&'static str, },
        { semantic_id:"Atom.isotope",item:callable,owner:type_,rust:crate::Atom::isotope,python:"isotope",javascript:"isotope",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<u16>,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->Option<u16>, },
        { semantic_id:"Atom.atom_map",item:callable,owner:type_,rust:crate::Atom::atom_map,python:"atom_map",javascript:"atomMap",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<u32>,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->Option<u32>, },
        { semantic_id:"Atom.is_aromatic",item:callable,owner:type_,rust:crate::Atom::is_aromatic,python:"is_aromatic",javascript:"isAromatic",feature:"runtime",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->bool, },
        { semantic_id:"Atom.explicit_hydrogens",item:callable,owner:type_,rust:crate::Atom::explicit_hydrogens,python:"explicit_hydrogens",javascript:"explicitHydrogens",feature:"runtime",status:experimental,kind:instance,parameters:[],output:u8,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->u8, },
        { semantic_id:"Atom.no_implicit",item:callable,owner:type_,rust:crate::Atom::no_implicit,python:"no_implicit",javascript:"noImplicit",feature:"runtime",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->bool, },
        { semantic_id:"Atom.radical_electrons",item:callable,owner:type_,rust:crate::Atom::radical_electrons,python:"radical_electrons",javascript:"radicalElectrons",feature:"runtime",status:experimental,kind:instance,parameters:[],output:u8,error:none,state:read_only,operation:none,signature:fn(&crate::Atom)->u8, },
        { semantic_id:"Bond.id",item:callable,owner:type_,rust:crate::Bond::id,python:"id",javascript:"id",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::BondId,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::BondId, },
        { semantic_id:"Bond.begin",item:callable,owner:type_,rust:crate::Bond::begin,python:"begin",javascript:"begin",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::AtomId,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::AtomId, },
        { semantic_id:"Bond.end",item:callable,owner:type_,rust:crate::Bond::end,python:"end",javascript:"end",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::AtomId,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::AtomId, },
        { semantic_id:"Bond.order",item:callable,owner:type_,rust:crate::Bond::order,python:"order",javascript:"order",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::BondOrder,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::BondOrder, },
        { semantic_id:"Bond.order_code",item:callable,owner:type_,rust:crate::Bond::order_code,python:"order_code",javascript:"orderCode",feature:"runtime",status:experimental,kind:instance,parameters:[],output:i64,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->i64, },
        { semantic_id:"Bond.order_name",item:callable,owner:type_,rust:crate::Bond::order_name,python:"order_name",javascript:"orderName",feature:"runtime",status:experimental,kind:instance,parameters:[],output:&'static str,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->&'static str, },
        { semantic_id:"Bond.direction",item:callable,owner:type_,rust:crate::Bond::direction,python:"direction",javascript:"direction",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::BondDirection,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::BondDirection, },
        { semantic_id:"Bond.direction_code",item:callable,owner:type_,rust:crate::Bond::direction_code,python:"direction_code",javascript:"directionCode",feature:"runtime",status:experimental,kind:instance,parameters:[],output:i64,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->i64, },
        { semantic_id:"Bond.direction_name",item:callable,owner:type_,rust:crate::Bond::direction_name,python:"direction_name",javascript:"directionName",feature:"runtime",status:experimental,kind:instance,parameters:[],output:&'static str,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->&'static str, },
        { semantic_id:"Bond.stereo",item:callable,owner:type_,rust:crate::Bond::stereo,python:"stereo",javascript:"stereo",feature:"runtime",status:experimental,kind:instance,parameters:[],output:crate::BondStereo,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->crate::BondStereo, },
        { semantic_id:"Bond.stereo_code",item:callable,owner:type_,rust:crate::Bond::stereo_code,python:"stereo_code",javascript:"stereoCode",feature:"runtime",status:experimental,kind:instance,parameters:[],output:i64,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->i64, },
        { semantic_id:"Bond.stereo_name",item:callable,owner:type_,rust:crate::Bond::stereo_name,python:"stereo_name",javascript:"stereoName",feature:"runtime",status:experimental,kind:instance,parameters:[],output:&'static str,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->&'static str, },
        { semantic_id:"Bond.stereo_atoms",item:callable,owner:type_,rust:crate::Bond::stereo_atoms,python:"stereo_atoms",javascript:"stereoAtoms",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<[crate::AtomId; 2]>,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->Option<[crate::AtomId; 2]>, },
        { semantic_id:"Bond.is_aromatic",item:callable,owner:type_,rust:crate::Bond::is_aromatic,python:"is_aromatic",javascript:"isAromatic",feature:"runtime",status:experimental,kind:instance,parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::Bond)->bool, },
        { semantic_id:"Atom.cip_descriptor",item:callable,owner:type_,rust:crate::Atom::cip_descriptor,python:"cip_descriptor",javascript:"cipDescriptor",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<crate::CipDescriptor>,error:crate::CipDescriptorError,state:read_only,operation:none,signature:fn(&crate::Atom)->Result<Option<crate::CipDescriptor>,crate::CipDescriptorError>, },
        { semantic_id:"Atom.cip_neighbor_order",item:callable,owner:type_,rust:crate::Atom::cip_neighbor_order,python:"cip_neighbor_order",javascript:"cipNeighborOrder",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<Vec<u32>>,error:crate::CipDescriptorError,state:read_only,operation:none,signature:fn(&crate::Atom)->Result<Option<Vec<u32>>,crate::CipDescriptorError>, },
        { semantic_id:"Atom.cip_rank",item:callable,owner:type_,rust:crate::Atom::cip_rank,python:"cip_rank",javascript:"cipRank",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<u32>,error:crate::PropertyValueError,state:read_only,operation:none,signature:fn(&crate::Atom)->Result<Option<u32>,crate::PropertyValueError>, },
        { semantic_id:"Bond.cip_descriptor",item:callable,owner:type_,rust:crate::Bond::cip_descriptor,python:"cip_descriptor",javascript:"cipDescriptor",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<crate::CipDescriptor>,error:crate::CipDescriptorError,state:read_only,operation:none,signature:fn(&crate::Bond)->Result<Option<crate::CipDescriptor>,crate::CipDescriptorError>, },
        { semantic_id:"Bond.cip_neighbor_order",item:callable,owner:type_,rust:crate::Bond::cip_neighbor_order,python:"cip_neighbor_order",javascript:"cipNeighborOrder",feature:"runtime",status:experimental,kind:instance,parameters:[],output:Option<Vec<u32>>,error:crate::CipDescriptorError,state:read_only,operation:none,signature:fn(&crate::Bond)->Result<Option<Vec<u32>>,crate::CipDescriptorError>, },
        #[cfg(feature="cap-stereo")]
        { semantic_id:"Molecule.cip_computed",item:callable,owner:molecule,rust:crate::Molecule::cip_computed,
          python:"cip_computed",javascript:"cipComputed",feature:"cap-stereo",status:experimental,kind:instance,
          parameters:[],output:bool,error:none,state:read_only,operation:none,signature:fn(&crate::Molecule)->bool, },
        {semantic_id:"AtomSpec.new",item:callable,owner:type_,rust:crate::AtomSpec::new,python:"new",javascript:"new",feature:"runtime",status:experimental,kind:static_,parameters:[{name:element,type:crate::Element,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::Element)->crate::AtomSpec,},
        {semantic_id:"BondSpec.new",item:callable,owner:type_,rust:crate::BondSpec::new,python:"new",javascript:"new",feature:"runtime",status:experimental,kind:static_,parameters:[{name:begin,type:crate::AtomId,default:required},{name:end,type:crate::AtomId,default:required},{name:order,type:crate::BondOrder,default:required},],output:crate::BondSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomId,crate::AtomId,crate::BondOrder)->crate::BondSpec,},
        {semantic_id:"AtomSpec.with_formal_charge",item:callable,owner:type_,rust:crate::AtomSpec::with_formal_charge,python:"with_formal_charge",javascript:"withFormalCharge",feature:"runtime",status:experimental,kind:instance,receiver:owned,parameters:[{name:value,type:i8,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomSpec,i8)->crate::AtomSpec,},
        {semantic_id:"AtomSpec.with_explicit_hydrogens",item:callable,owner:type_,rust:crate::AtomSpec::with_explicit_hydrogens,python:"with_explicit_hydrogens",javascript:"withExplicitHydrogens",feature:"runtime",status:experimental,kind:instance,receiver:owned,parameters:[{name:value,type:u8,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomSpec,u8)->crate::AtomSpec,},
        {semantic_id:"AtomSpec.with_atom_map",item:callable,owner:type_,rust:crate::AtomSpec::with_atom_map,python:"with_atom_map",javascript:"withAtomMap",feature:"runtime",status:experimental,kind:instance,receiver:owned,parameters:[{name:value,type:u32,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomSpec,u32)->crate::AtomSpec,},
        {semantic_id:"AtomSpec.with_isotope",item:callable,owner:type_,rust:crate::AtomSpec::with_isotope,python:"with_isotope",javascript:"withIsotope",feature:"runtime",status:experimental,kind:instance,receiver:owned,parameters:[{name:value,type:u16,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomSpec,u16)->crate::AtomSpec,},
        {semantic_id:"AtomSpec.with_no_implicit",item:callable,owner:type_,rust:crate::AtomSpec::with_no_implicit,python:"with_no_implicit",javascript:"withNoImplicit",feature:"runtime",status:experimental,kind:instance,receiver:owned,parameters:[{name:value,type:bool,default:required},],output:crate::AtomSpec,error:none,state:value_returning,operation:none,signature:fn(crate::AtomSpec,bool)->crate::AtomSpec,},
        #[cfg(feature="cap-sanitize")]
        {semantic_id:"Molecule.sanitize_",item:callable,owner:molecule,rust:crate::Molecule::sanitize_,python:"sanitize_",javascript:"sanitize_",feature:"cap-sanitize",kind:instance,parameters:[],output:(),error:crate::OperationError,state:in_place,operation:"sanitize_",signature:fn(&mut crate::Molecule)->Result<(),crate::OperationError>,},
        #[cfg(feature="cap-sanitize")]
        {semantic_id:"Molecule.sanitize_with_params_",item:callable,owner:molecule,rust:crate::Molecule::sanitize_with_params_,python:"sanitize_with_params_",javascript:"sanitizeWithParams_",feature:"cap-sanitize",kind:instance,parameters:[{name:params,type:&crate::SanitizeParams,default:required},],output:(),error:crate::OperationError,state:in_place,operation:"sanitize_with_params_",signature:for<'a,'b> fn(&'a mut crate::Molecule,&'b crate::SanitizeParams)->Result<(),crate::OperationError>,},
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "types.LegacyTopologicalTorsionParams", item: type, owner: type_, rust: crate::LegacyTopologicalTorsionParams, python: "LegacyTopologicalTorsionParams", javascript: "LegacyTopologicalTorsionParams", feature: "cap-fingerprints", status: experimental, role: parameter, },
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.TopologicalTorsionFingerprintGenerator",item:type,owner:type_,rust:crate::TopologicalTorsionFingerprintGenerator,python:"TopologicalTorsionFingerprintGenerator",javascript:"TopologicalTorsionFingerprintGenerator",feature:"cap-fingerprints",status:experimental,role:value},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.TopologicalTorsionSettings",item:type,owner:type_,rust:crate::TopologicalTorsionSettings,python:"TopologicalTorsionSettings",javascript:"TopologicalTorsionSettings",feature:"cap-fingerprints",status:experimental,role:value},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.TopologicalTorsionCallParams",item:type,owner:type_,rust:crate::TopologicalTorsionCallParams,python:"TopologicalTorsionCallParams",javascript:"TopologicalTorsionCallParams",feature:"cap-fingerprints",status:experimental,role:parameter},
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.new", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::new, python: "new", javascript: "new",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [{ name: params, type: Option<&crate::TopologicalTorsionParams>, default: none }, { name: atom_invariants, type: Option<crate::AtomPairAtomInvariantsGenerator>, default: none }], output: crate::TopologicalTorsionFingerprintGenerator, error: crate::TopologicalTorsionReadError,
            state: value_returning, operation: none, signature: fn(Option<&crate::TopologicalTorsionParams>,Option<crate::AtomPairAtomInvariantsGenerator>)->Result<crate::TopologicalTorsionFingerprintGenerator,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.from_json", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::from_json, python: "from_json", javascript: "fromJson",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [{ name: json, type: &str, default: required }], output: crate::TopologicalTorsionFingerprintGenerator, error: crate::TopologicalTorsionReadError,
            state: value_returning, operation: none, signature: fn(&str)->Result<crate::TopologicalTorsionFingerprintGenerator,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.settings", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::settings, python: "settings", javascript: "settings",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: crate::TopologicalTorsionSettings, error: none,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionFingerprintGenerator)->crate::TopologicalTorsionSettings,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.info_string", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::info_string, python: "info_string", javascript: "infoString",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: String, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionFingerprintGenerator)->Result<String,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.to_json", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::to_json, python: "to_json", javascript: "toJson",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: String, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionFingerprintGenerator)->Result<String,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.fingerprints", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::fingerprints, python: "fingerprints", javascript: "fingerprints",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: molecules, type: &[Option<&crate::Molecule>], default: required }, { name: num_threads, type: i32, default: literal(1) }], output: Vec<Option<crate::Fingerprint>>, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::TopologicalTorsionFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::Fingerprint>>,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.sparse_fingerprints", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::sparse_fingerprints, python: "sparse_fingerprints", javascript: "sparseFingerprints",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: molecules, type: &[Option<&crate::Molecule>], default: required }, { name: num_threads, type: i32, default: literal(1) }], output: Vec<Option<crate::SparseBitFingerprint>>, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::TopologicalTorsionFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseBitFingerprint>>,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.counts", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::counts, python: "counts", javascript: "counts",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: molecules, type: &[Option<&crate::Molecule>], default: required }, { name: num_threads, type: i32, default: literal(1) }], output: Vec<Option<crate::SparseCountFingerprint32>>, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::TopologicalTorsionFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseCountFingerprint32>>,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionFingerprintGenerator.sparse_counts", item: callable, owner: type_,
            rust: crate::TopologicalTorsionFingerprintGenerator::sparse_counts, python: "sparse_counts", javascript: "sparseCounts",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: molecules, type: &[Option<&crate::Molecule>], default: required }, { name: num_threads, type: i32, default: literal(1) }], output: Vec<Option<crate::SparseCountFingerprint>>, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::TopologicalTorsionFingerprintGenerator,&'b [Option<&'b crate::Molecule>],i32)->Result<Vec<Option<crate::SparseCountFingerprint>>,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.torsion_atom_count", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::torsion_atom_count, python: "torsion_atom_count", javascript: "torsionAtomCount",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: u32, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<u32,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_torsion_atom_count", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_torsion_atom_count, python: "set_torsion_atom_count", javascript: "setTorsionAtomCount",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: u32, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,u32)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.only_shortest_paths", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::only_shortest_paths, python: "only_shortest_paths", javascript: "onlyShortestPaths",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: bool, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<bool,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_only_shortest_paths", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_only_shortest_paths, python: "set_only_shortest_paths", javascript: "setOnlyShortestPaths",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: bool, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,bool)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.include_chirality", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::include_chirality, python: "include_chirality", javascript: "includeChirality",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: bool, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<bool,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_include_chirality", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_include_chirality, python: "set_include_chirality", javascript: "setIncludeChirality",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: bool, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,bool)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.count_simulation", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::count_simulation, python: "count_simulation", javascript: "countSimulation",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: bool, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<bool,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_count_simulation", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_count_simulation, python: "set_count_simulation", javascript: "setCountSimulation",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: bool, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,bool)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.fp_size", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::fp_size, python: "fp_size", javascript: "fpSize",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: u32, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<u32,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_fp_size", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_fp_size, python: "set_fp_size", javascript: "setFpSize",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: u32, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,u32)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.bits_per_feature", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::bits_per_feature, python: "bits_per_feature", javascript: "bitsPerFeature",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: u32, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<u32,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_bits_per_feature", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_bits_per_feature, python: "set_bits_per_feature", javascript: "setBitsPerFeature",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: u32, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,u32)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.count_bounds", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::count_bounds, python: "count_bounds", javascript: "countBounds",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: Vec<u32>, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<Vec<u32>,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.set_count_bounds", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::set_count_bounds, python: "set_count_bounds", javascript: "setCountBounds",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: mutable,
            parameters: [{ name: value, type: Vec<u32>, default: required }], output: (), error: crate::TopologicalTorsionReadError,
            state: in_place, operation: none, signature: for<'a> fn(&'a mut crate::TopologicalTorsionSettings,Vec<u32>)->Result<(),crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "TopologicalTorsionSettings.params", item: callable, owner: type_,
            rust: crate::TopologicalTorsionSettings::params, python: "params", javascript: "params",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [], output: crate::TopologicalTorsionParams, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a> fn(&'a crate::TopologicalTorsionSettings)->Result<crate::TopologicalTorsionParams,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.topological_torsion_fingerprint_with_generator", item: callable, owner: molecule,
            rust: crate::Molecule::topological_torsion_fingerprint_with_generator, python: "topological_torsion_fingerprint_with_generator", javascript: "topologicalTorsionFingerprintWithGenerator",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: generator, type: &crate::TopologicalTorsionFingerprintGenerator, default: required }, { name: params, type: Option<&crate::TopologicalTorsionCallParams>, default: none }, { name: output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: none }], output: crate::Fingerprint, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::TopologicalTorsionFingerprintGenerator,Option<&'c crate::TopologicalTorsionCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::Fingerprint,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.topological_torsion_sparse_fingerprint_with_generator", item: callable, owner: molecule,
            rust: crate::Molecule::topological_torsion_sparse_fingerprint_with_generator, python: "topological_torsion_sparse_fingerprint_with_generator", javascript: "topologicalTorsionSparseFingerprintWithGenerator",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: generator, type: &crate::TopologicalTorsionFingerprintGenerator, default: required }, { name: params, type: Option<&crate::TopologicalTorsionCallParams>, default: none }, { name: output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: none }], output: crate::SparseBitFingerprint, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::TopologicalTorsionFingerprintGenerator,Option<&'c crate::TopologicalTorsionCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseBitFingerprint,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.topological_torsion_count_fingerprint_with_generator", item: callable, owner: molecule,
            rust: crate::Molecule::topological_torsion_count_fingerprint_with_generator, python: "topological_torsion_count_fingerprint_with_generator", javascript: "topologicalTorsionCountFingerprintWithGenerator",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: generator, type: &crate::TopologicalTorsionFingerprintGenerator, default: required }, { name: params, type: Option<&crate::TopologicalTorsionCallParams>, default: none }, { name: output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: none }], output: crate::SparseCountFingerprint32, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::TopologicalTorsionFingerprintGenerator,Option<&'c crate::TopologicalTorsionCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseCountFingerprint32,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Molecule.topological_torsion_sparse_count_fingerprint_with_generator", item: callable, owner: molecule,
            rust: crate::Molecule::topological_torsion_sparse_count_fingerprint_with_generator, python: "topological_torsion_sparse_count_fingerprint_with_generator", javascript: "topologicalTorsionSparseCountFingerprintWithGenerator",
            feature: "cap-fingerprints", status: experimental, kind: instance,
            parameters: [{ name: generator, type: &crate::TopologicalTorsionFingerprintGenerator, default: required }, { name: params, type: Option<&crate::TopologicalTorsionCallParams>, default: none }, { name: output, type: Option<&mut crate::FingerprintAdditionalOutput>, default: none }], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError,
            state: read_only, operation: none, signature: for<'a,'b,'c,'d> fn(&'a crate::Molecule,&'b crate::TopologicalTorsionFingerprintGenerator,Option<&'c crate::TopologicalTorsionCallParams>,Option<&'d mut crate::FingerprintAdditionalOutput>)->Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "LegacyTopologicalTorsionParams.new", item: callable, owner: type_, rust: crate::LegacyTopologicalTorsionParams::new, python: "new", javascript: "new", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [{ name: torsion_atom_count, type: u32, default: integer(4) },{ name: include_chirality, type: bool, default: boolean(false) },{ name: fp_size, type: u32, default: integer(2048) },{ name: bits_per_entry, type: u32, default: integer(4) },{ name: from_atoms, type: Option<Vec<u32>>, default: none },{ name: ignore_atoms, type: Option<Vec<u32>>, default: none },{ name: custom_atom_invariants, type: Option<Vec<u32>>, default: none }], output: crate::LegacyTopologicalTorsionParams, error: none, state: value_returning, operation: none, signature: fn(u32,bool,u32,u32,Option<Vec<u32>>,Option<Vec<u32>>,Option<Vec<u32>>)->crate::LegacyTopologicalTorsionParams, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_sparse_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_sparse_count_fingerprint, python: "legacy_topological_torsion_sparse_count_fingerprint", javascript: "legacyTopologicalTorsionSparseCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule)->Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_sparse_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_sparse_count_fingerprint_with_params, python: "legacy_topological_torsion_sparse_count_fingerprint_with_params", javascript: "legacyTopologicalTorsionSparseCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::LegacyTopologicalTorsionParams, default: required }], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::Molecule,&'b crate::LegacyTopologicalTorsionParams)->Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_count_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_count_fingerprint, python: "legacy_topological_torsion_count_fingerprint", javascript: "legacyTopologicalTorsionCountFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule)->Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_count_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_count_fingerprint_with_params, python: "legacy_topological_torsion_count_fingerprint_with_params", javascript: "legacyTopologicalTorsionCountFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::LegacyTopologicalTorsionParams, default: required }], output: crate::SparseCountFingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::Molecule,&'b crate::LegacyTopologicalTorsionParams)->Result<crate::SparseCountFingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_fingerprint", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_fingerprint, python: "legacy_topological_torsion_fingerprint", javascript: "legacyTopologicalTorsionFingerprint", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: crate::Fingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule)->Result<crate::Fingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.legacy_topological_torsion_fingerprint_with_params", item: callable, owner: molecule, rust: crate::Molecule::legacy_topological_torsion_fingerprint_with_params, python: "legacy_topological_torsion_fingerprint_with_params", javascript: "legacyTopologicalTorsionFingerprintWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: params, type: &crate::LegacyTopologicalTorsionParams, default: required }], output: crate::Fingerprint, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: for<'a,'b> fn(&'a crate::Molecule,&'b crate::LegacyTopologicalTorsionParams)->Result<crate::Fingerprint,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_ids", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_ids, python: "topological_torsion_ids", javascript: "topologicalTorsionIds", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: Vec<u64>, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule)->Result<Vec<u64>,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "Molecule.topological_torsion_ids_with_params", item: callable, owner: molecule, rust: crate::Molecule::topological_torsion_ids_with_params, python: "topological_torsion_ids_with_params", javascript: "topologicalTorsionIdsWithParams", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [{ name: torsion_atom_count, type: u32, default: required }], output: Vec<u64>, error: crate::TopologicalTorsionReadError, state: read_only, operation: none, signature: fn(&crate::Molecule,u32)->Result<Vec<u64>,crate::TopologicalTorsionReadError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomCodeExplanation.from_code", item: callable, owner: type_, rust: crate::AtomCodeExplanation::from_code, python: "from_code", javascript: "fromCode", feature: "cap-fingerprints", status: experimental, kind: static_, parameters: [{ name: code, type: u64, default: required }, { name: branch_subtract, type: i64, default: integer(0) }, { name: include_chirality, type: bool, default: boolean(false) }], output: crate::AtomCodeExplanation, error: crate::AtomCodeExplanationError, state: value_returning, operation: none, signature: fn(u64,i64,bool)->Result<crate::AtomCodeExplanation,crate::AtomCodeExplanationError>, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomCodeExplanation.symbol", item: callable, owner: type_, rust: crate::AtomCodeExplanation::symbol, python: "symbol", javascript: "symbol", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: &'static str, error: none, state: read_only, operation: none, signature: fn(&crate::AtomCodeExplanation)->&'static str, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomCodeExplanation.branch_count", item: callable, owner: type_, rust: crate::AtomCodeExplanation::branch_count, python: "branch_count", javascript: "branchCount", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn(&crate::AtomCodeExplanation)->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomCodeExplanation.pi_electrons", item: callable, owner: type_, rust: crate::AtomCodeExplanation::pi_electrons, python: "pi_electrons", javascript: "piElectrons", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: u32, error: none, state: read_only, operation: none, signature: fn(&crate::AtomCodeExplanation)->u32, },
        #[cfg(feature = "cap-fingerprints")]
        { semantic_id: "AtomCodeExplanation.chirality", item: callable, owner: type_, rust: crate::AtomCodeExplanation::chirality, python: "chirality", javascript: "chirality", feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared, parameters: [], output: Option<&'static str>, error: none, state: read_only, operation: none, signature: fn(&crate::AtomCodeExplanation)->Option<&'static str>, },
        #[cfg(feature = "cap-serialization")]
        { semantic_id:"types.PickleError", item:type, owner:type_,
          rust:crate::PickleError, python:"PickleError", javascript:"PickleError",
          feature:"cap-serialization", status:experimental, role:error, },
        #[cfg(feature = "cap-serialization")]
        { semantic_id:"Molecule.to_binary", item:callable, owner:molecule,
          rust:crate::Molecule::to_binary, python:"to_binary", javascript:"toBinary",
          feature:"cap-serialization", status:experimental, kind:instance,
          parameters:[], output:Vec<u8>, error:crate::PickleError,
          state:read_only, operation:none,
          signature:fn(&crate::Molecule)->Result<Vec<u8>,crate::PickleError>, },
        #[cfg(feature = "cap-serialization")]
        { semantic_id:"Molecule.from_binary", item:callable, owner:molecule,
          rust:crate::Molecule::from_binary, python:"from_binary", javascript:"fromBinary",
          feature:"cap-serialization", status:experimental, kind:static_,
          parameters:[{name:data,type:&[u8],default:required}], output:crate::Molecule, error:crate::PickleError,
          state:value_returning, operation:none,
          signature:fn(&[u8])->Result<crate::Molecule,crate::PickleError>, },
        #[cfg(feature = "cap-hashing")]
        { semantic_id:"types.CipRankError", item:type, owner:type_, rust:crate::CipRankError, python:"CipRankError", javascript:"CipRankError", feature:"cap-hashing", status:experimental, role:error, },
        #[cfg(feature = "cap-hashing")]
        { semantic_id:"types.MoleculeHashError", item:type, owner:type_, rust:crate::MoleculeHashError,
          python:"MoleculeHashError", javascript:"MoleculeHashError", feature:"cap-hashing", status:experimental, role:error, },
        #[cfg(feature = "cap-hashing")]
        { semantic_id:"Molecule.molecular_hash", item:callable, owner:molecule, rust:crate::Molecule::molecular_hash,
          python:"molecular_hash", javascript:"molecularHash", feature:"cap-hashing", status:experimental, kind:instance,
          parameters:[], output:u64, error:crate::MoleculeHashError, state:read_only, operation:none,
          signature:fn(&crate::Molecule)->Result<u64,crate::MoleculeHashError>, },
        #[cfg(feature = "cap-hashing")]
        { semantic_id:"Molecule.molecular_hash_with_ranks", item:callable, owner:molecule, rust:crate::Molecule::molecular_hash_with_ranks,
          python:"molecular_hash_with_ranks", javascript:"molecularHashWithRanks", feature:"cap-hashing", status:experimental, kind:instance,
          parameters:[{name:ranks,type:&[u32],default:required}], output:u64, error:crate::MoleculeHashError, state:read_only, operation:none,
          signature:fn(&crate::Molecule,&[u32])->Result<u64,crate::MoleculeHashError>, },

        { semantic_id: "types.LigandRef", item: type, owner: type_, rust: crate::LigandRef, python: "LigandRef", javascript: "LigandRef", feature: "runtime", status: native, role: value, },
        { semantic_id: "types.TetrahedralStereo", item: type, owner: type_, rust: crate::TetrahedralStereo, python: "TetrahedralStereo", javascript: "TetrahedralStereo", feature: "runtime", status: native, role: result, },
        #[cfg(feature = "cap-stereo")]
        { semantic_id: "types.StereoReadError", item: type, owner: type_, rust: crate::StereoReadError, python: "StereoReadError", javascript: "StereoReadError", feature: "cap-stereo", status: native, role: error, },
        #[cfg(feature = "cap-stereo")]
        { semantic_id: "Molecule.tetrahedral_stereo", item: callable, owner: molecule, rust: crate::Molecule::tetrahedral_stereo, python: "tetrahedral_stereo", javascript: "tetrahedralStereo", feature: "cap-stereo", status: native, kind: instance, parameters: [], output: Vec<crate::TetrahedralStereo>, error: crate::StereoReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<Vec<crate::TetrahedralStereo>, crate::StereoReadError>, },
        #[cfg(feature = "cap-stereo")]
        { semantic_id: "Molecule.perceive_stereochemistry", item: callable, owner: molecule, rust: crate::Molecule::perceive_stereochemistry, python: "perceive_stereochemistry", javascript: "perceiveStereochemistry", feature: "cap-stereo", status: native, kind: instance, parameters: [], output: (), error: crate::StereoReadError, state: read_only, operation: none, signature: fn(&crate::Molecule) -> Result<(), crate::StereoReadError>, },
        #[cfg(feature = "cap-stereo")]
        { semantic_id: "Molecule.find_chiral_centers", item: callable, owner: molecule, rust: crate::Molecule::find_chiral_centers, python: "find_chiral_centers", javascript: "findChiralCenters", feature: "cap-stereo", status: native, kind: instance, parameters: [{name: include_unassigned, type: bool, default: true}], output: Vec<(usize, String)>, error: none, state: read_only, operation: none, signature: fn(&crate::Molecule, bool) -> Vec<(usize, String)>, },





#[cfg(feature="cap-alignment")]
{semantic_id:"types.AlignmentAtomMap",item:type,owner:type_,rust:crate::AlignmentAtomMap,python:"AlignmentAtomMap",javascript:"AlignmentAtomMap",feature:"cap-alignment",status:experimental,role:value,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.AlignmentParameters",item:type,owner:type_,rust:crate::AlignmentParameters,python:"AlignmentParameters",javascript:"AlignmentParameters",feature:"cap-alignment",status:experimental,role:parameter,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.BestAlignmentParameters",item:type,owner:type_,rust:crate::BestAlignmentParameters,python:"BestAlignmentParameters",javascript:"BestAlignmentParameters",feature:"cap-alignment",status:experimental,role:parameter,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.CoordinateRmsdParameters",item:type,owner:type_,rust:crate::CoordinateRmsdParameters,python:"CoordinateRmsdParameters",javascript:"CoordinateRmsdParameters",feature:"cap-alignment",status:experimental,role:parameter,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.AllConformerRmsdParameters",item:type,owner:type_,rust:crate::AllConformerRmsdParameters,python:"AllConformerRmsdParameters",javascript:"AllConformerRmsdParameters",feature:"cap-alignment",status:experimental,role:parameter,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.ConformerAlignmentParameters",item:type,owner:type_,rust:crate::ConformerAlignmentParameters,python:"ConformerAlignmentParameters",javascript:"ConformerAlignmentParameters",feature:"cap-alignment",status:experimental,role:parameter,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.AlignmentResult",item:type,owner:type_,rust:crate::AlignmentResult,python:"AlignmentResult",javascript:"AlignmentResult",feature:"cap-alignment",status:experimental,role:result,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.AlignmentTransform",item:type,owner:type_,rust:crate::AlignmentTransform,python:"AlignmentTransform",javascript:"AlignmentTransform",feature:"cap-alignment",status:experimental,role:result,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.ConformerRmsd",item:type,owner:type_,rust:crate::ConformerRmsd,python:"ConformerRmsd",javascript:"ConformerRmsd",feature:"cap-alignment",status:experimental,role:result,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.ConformerAlignmentReport",item:type,owner:type_,rust:crate::ConformerAlignmentReport,python:"ConformerAlignmentReport",javascript:"ConformerAlignmentReport",feature:"cap-alignment",status:experimental,role:result,},
#[cfg(feature="cap-alignment")]
{semantic_id:"types.AlignmentError",item:type,owner:type_,rust:crate::AlignmentError,python:"AlignmentError",javascript:"AlignmentError",feature:"cap-alignment",status:experimental,role:error,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.alignment_transform_to",item:callable,owner:molecule,rust:crate::Molecule::alignment_transform_to,python:"alignment_transform_to",javascript:"alignmentTransformTo",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:crate::AlignmentResult,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule)->Result<crate::AlignmentResult,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.alignment_transform_to_with_params",item:callable,owner:molecule,rust:crate::Molecule::alignment_transform_to_with_params,python:"alignment_transform_to_with_params",javascript:"alignmentTransformToWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::AlignmentParameters,default:required}],output:crate::AlignmentResult,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule,&crate::AlignmentParameters)->Result<crate::AlignmentResult,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.best_alignment_to",item:callable,owner:molecule,rust:crate::Molecule::best_alignment_to,python:"best_alignment_to",javascript:"bestAlignmentTo",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:crate::AlignmentResult,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule)->Result<crate::AlignmentResult,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.best_alignment_to_with_params",item:callable,owner:molecule,rust:crate::Molecule::best_alignment_to_with_params,python:"best_alignment_to_with_params",javascript:"bestAlignmentToWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::BestAlignmentParameters,default:required}],output:crate::AlignmentResult,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule,&crate::BestAlignmentParameters)->Result<crate::AlignmentResult,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.best_rmsd_to",item:callable,owner:molecule,rust:crate::Molecule::best_rmsd_to,python:"best_rmsd_to",javascript:"bestRmsdTo",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:f64,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule)->Result<f64,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.best_rmsd_to_with_params",item:callable,owner:molecule,rust:crate::Molecule::best_rmsd_to_with_params,python:"best_rmsd_to_with_params",javascript:"bestRmsdToWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::BestAlignmentParameters,default:required}],output:f64,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule,&crate::BestAlignmentParameters)->Result<f64,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.coordinate_rmsd_to",item:callable,owner:molecule,rust:crate::Molecule::coordinate_rmsd_to,python:"coordinate_rmsd_to",javascript:"coordinateRmsdTo",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:f64,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule)->Result<f64,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.coordinate_rmsd_to_with_params",item:callable,owner:molecule,rust:crate::Molecule::coordinate_rmsd_to_with_params,python:"coordinate_rmsd_to_with_params",javascript:"coordinateRmsdToWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::CoordinateRmsdParameters,default:required}],output:f64,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::Molecule,&crate::CoordinateRmsdParameters)->Result<f64,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.all_conformer_best_rmsds",item:callable,owner:molecule,rust:crate::Molecule::all_conformer_best_rmsds,python:"all_conformer_best_rmsds",javascript:"allConformerBestRmsds",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[],output:Vec<crate::ConformerRmsd>,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule)->Result<Vec<crate::ConformerRmsd>,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.all_conformer_best_rmsds_with_params",item:callable,owner:molecule,rust:crate::Molecule::all_conformer_best_rmsds_with_params,python:"all_conformer_best_rmsds_with_params",javascript:"allConformerBestRmsdsWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:params,type:&crate::AllConformerRmsdParameters,default:required}],output:Vec<crate::ConformerRmsd>,error:crate::AlignmentError,state:read_only,operation:none,signature:fn(&crate::Molecule,&crate::AllConformerRmsdParameters)->Result<Vec<crate::ConformerRmsd>,crate::AlignmentError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.with_alignment_to",item:callable,owner:molecule,rust:crate::Molecule::with_alignment_to,python:"with_alignment_to",javascript:"withAlignmentTo",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:(crate::Molecule,crate::AlignmentResult),error:crate::OperationError,state:value_returning,operation:"with_alignment_to",signature:fn(&crate::Molecule,&crate::Molecule)->Result<(crate::Molecule,crate::AlignmentResult),crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.with_alignment_to_with_params",item:callable,owner:molecule,rust:crate::Molecule::with_alignment_to_with_params,python:"with_alignment_to_with_params",javascript:"withAlignmentToWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::AlignmentParameters,default:required}],output:(crate::Molecule,crate::AlignmentResult),error:crate::OperationError,state:value_returning,operation:"with_alignment_to_with_params",signature:fn(&crate::Molecule,&crate::Molecule,&crate::AlignmentParameters)->Result<(crate::Molecule,crate::AlignmentResult),crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.align_to_",item:callable,owner:molecule,rust:crate::Molecule::align_to_,python:"align_to_",javascript:"alignTo",feature:"cap-alignment",kind:instance,receiver:mutable,parameters:[{name:reference,type:&crate::Molecule,default:required}],output:crate::AlignmentResult,error:crate::OperationError,state:in_place,operation:"align_to_",signature:fn(&mut crate::Molecule,&crate::Molecule)->Result<crate::AlignmentResult,crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.align_to_with_params_",item:callable,owner:molecule,rust:crate::Molecule::align_to_with_params_,python:"align_to_with_params_",javascript:"alignToWithParams",feature:"cap-alignment",kind:instance,receiver:mutable,parameters:[{name:reference,type:&crate::Molecule,default:required},{name:params,type:&crate::AlignmentParameters,default:required}],output:crate::AlignmentResult,error:crate::OperationError,state:in_place,operation:"align_to_with_params_",signature:fn(&mut crate::Molecule,&crate::Molecule,&crate::AlignmentParameters)->Result<crate::AlignmentResult,crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.with_aligned_conformers",item:callable,owner:molecule,rust:crate::Molecule::with_aligned_conformers,python:"with_aligned_conformers",javascript:"withAlignedConformers",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[],output:(crate::Molecule,crate::ConformerAlignmentReport),error:crate::OperationError,state:value_returning,operation:"with_aligned_conformers",signature:fn(&crate::Molecule)->Result<(crate::Molecule,crate::ConformerAlignmentReport),crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.with_aligned_conformers_with_params",item:callable,owner:molecule,rust:crate::Molecule::with_aligned_conformers_with_params,python:"with_aligned_conformers_with_params",javascript:"withAlignedConformersWithParams",feature:"cap-alignment",kind:instance,receiver:shared,parameters:[{name:params,type:&crate::ConformerAlignmentParameters,default:required}],output:(crate::Molecule,crate::ConformerAlignmentReport),error:crate::OperationError,state:value_returning,operation:"with_aligned_conformers_with_params",signature:fn(&crate::Molecule,&crate::ConformerAlignmentParameters)->Result<(crate::Molecule,crate::ConformerAlignmentReport),crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.align_conformers_",item:callable,owner:molecule,rust:crate::Molecule::align_conformers_,python:"align_conformers_",javascript:"alignConformers",feature:"cap-alignment",kind:instance,receiver:mutable,parameters:[],output:crate::ConformerAlignmentReport,error:crate::OperationError,state:in_place,operation:"align_conformers_",signature:fn(&mut crate::Molecule)->Result<crate::ConformerAlignmentReport,crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"Molecule.align_conformers_with_params_",item:callable,owner:molecule,rust:crate::Molecule::align_conformers_with_params_,python:"align_conformers_with_params_",javascript:"alignConformersWithParams",feature:"cap-alignment",kind:instance,receiver:mutable,parameters:[{name:params,type:&crate::ConformerAlignmentParameters,default:required}],output:crate::ConformerAlignmentReport,error:crate::OperationError,state:in_place,operation:"align_conformers_with_params_",signature:fn(&mut crate::Molecule,&crate::ConformerAlignmentParameters)->Result<crate::ConformerAlignmentReport,crate::OperationError>,},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentResult.rmsd",item:callable,owner:type_,rust:crate::AlignmentResult::rmsd,python:"rmsd",javascript:"rmsd",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::AlignmentResult)->f64,},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentResult.transform",item:callable,owner:type_,rust:crate::AlignmentResult::transform,python:"transform",javascript:"transform",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:&'a crate::AlignmentTransform,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::AlignmentResult)->&'a crate::AlignmentTransform,},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentResult.atom_map",item:callable,owner:type_,rust:crate::AlignmentResult::atom_map,python:"atom_map",javascript:"atomMap",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:&'a [crate::AlignmentAtomMap],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::AlignmentResult)->&'a [crate::AlignmentAtomMap],},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentTransform.matrix",item:callable,owner:type_,rust:crate::AlignmentTransform::matrix,python:"matrix",javascript:"matrix",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:&'a [[f64;4];4],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::AlignmentTransform)->&'a [[f64;4];4],},
#[cfg(feature="cap-alignment")]
{semantic_id:"ConformerRmsd.rmsd",item:callable,owner:type_,rust:crate::ConformerRmsd::rmsd,python:"rmsd",javascript:"rmsd",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:f64,error:none,state:read_only,operation:none,signature:fn(&crate::ConformerRmsd)->f64,},
#[cfg(feature="cap-alignment")]
{semantic_id:"ConformerRmsd.probe_conformer_id",item:callable,owner:type_,rust:crate::ConformerRmsd::probe_conformer_id,python:"probe_conformer_id",javascript:"probeConformerId",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:usize,error:none,state:read_only,operation:none,signature:fn(&crate::ConformerRmsd)->usize,},
#[cfg(feature="cap-alignment")]
{semantic_id:"ConformerRmsd.reference_conformer_id",item:callable,owner:type_,rust:crate::ConformerRmsd::reference_conformer_id,python:"reference_conformer_id",javascript:"referenceConformerId",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:usize,error:none,state:read_only,operation:none,signature:fn(&crate::ConformerRmsd)->usize,},
#[cfg(feature="cap-alignment")]
{semantic_id:"ConformerAlignmentReport.rmsds",item:callable,owner:type_,rust:crate::ConformerAlignmentReport::rmsds,python:"rmsds",javascript:"rmsds",feature:"cap-alignment",status:experimental,kind:instance,parameters:[],output:&'a [f64],error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::ConformerAlignmentReport)->&'a [f64],},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentAtomMap.new",item:callable,owner:type_,rust:crate::AlignmentAtomMap::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:probe_atom,type: usize,default:required},{name:reference_atom,type: usize,default:required}],output:crate::AlignmentAtomMap,error:none,state:value_returning,operation:none,signature:fn( usize, usize)->crate::AlignmentAtomMap,},
#[cfg(feature="cap-alignment")]
{semantic_id:"AlignmentParameters.new",item:callable,owner:type_,rust:crate::AlignmentParameters::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:probe_conformer_id,type:i32,default:integer(-1)},{name:reference_conformer_id,type:i32,default:integer(-1)},{name:atom_map,type:Option<Vec<crate::AlignmentAtomMap>>,default:none},{name:weights,type:Option<Vec<f64>>,default:none},{name:reflect,type:bool,default:boolean(false)},{name:max_iterations,type:u32,default:integer(50)}],output:crate::AlignmentParameters,error:none,state:value_returning,operation:none,signature:fn(i32,i32,Option<Vec<crate::AlignmentAtomMap>>,Option<Vec<f64>>,bool,u32)->crate::AlignmentParameters,},
#[cfg(feature="cap-alignment")]
{semantic_id:"BestAlignmentParameters.new",item:callable,owner:type_,rust:crate::BestAlignmentParameters::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:probe_conformer_id,type:i32,default:integer(-1)},{name:reference_conformer_id,type:i32,default:integer(-1)},{name:atom_maps,type:Vec<Vec<crate::AlignmentAtomMap>>,default:"[]"},{name:weights,type:Option<Vec<f64>>,default:none},{name:reflect,type:bool,default:boolean(false)},{name:max_iterations,type:u32,default:integer(50)},{name:max_matches,type:i32,default:integer(1000000)},{name:symmetrize_conjugated_terminal_groups,type:bool,default:boolean(true)},{name:ignore_hydrogens,type:bool,default:boolean(true)},{name:num_threads,type:i32,default:integer(1)}],output:crate::BestAlignmentParameters,error:none,state:value_returning,operation:none,signature:fn(i32,i32,Vec<Vec<crate::AlignmentAtomMap>>,Option<Vec<f64>>,bool,u32,i32,bool,bool,i32)->crate::BestAlignmentParameters,},
#[cfg(feature="cap-alignment")]
{semantic_id:"CoordinateRmsdParameters.new",item:callable,owner:type_,rust:crate::CoordinateRmsdParameters::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:probe_conformer_id,type:i32,default:integer(-1)},{name:reference_conformer_id,type:i32,default:integer(-1)},{name:atom_maps,type:Vec<Vec<crate::AlignmentAtomMap>>,default:"[]"},{name:weights,type:Option<Vec<f64>>,default:none},{name:max_matches,type:i32,default:integer(1000000)},{name:symmetrize_conjugated_terminal_groups,type:bool,default:boolean(true)}],output:crate::CoordinateRmsdParameters,error:none,state:value_returning,operation:none,signature:fn(i32,i32,Vec<Vec<crate::AlignmentAtomMap>>,Option<Vec<f64>>,i32,bool)->crate::CoordinateRmsdParameters,},
#[cfg(feature="cap-alignment")]
{semantic_id:"AllConformerRmsdParameters.new",item:callable,owner:type_,rust:crate::AllConformerRmsdParameters::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:atom_maps,type:Vec<Vec<crate::AlignmentAtomMap>>,default:"[]"},{name:weights,type:Option<Vec<f64>>,default:none},{name:max_matches,type:i32,default:integer(1000000)},{name:symmetrize_conjugated_terminal_groups,type:bool,default:boolean(true)},{name:ignore_hydrogens,type:bool,default:boolean(true)},{name:num_threads,type:i32,default:integer(1)}],output:crate::AllConformerRmsdParameters,error:none,state:value_returning,operation:none,signature:fn(Vec<Vec<crate::AlignmentAtomMap>>,Option<Vec<f64>>,i32,bool,bool,i32)->crate::AllConformerRmsdParameters,},
#[cfg(feature="cap-alignment")]
{semantic_id:"ConformerAlignmentParameters.new",item:callable,owner:type_,rust:crate::ConformerAlignmentParameters::new,python:"new",javascript:"new",feature:"cap-alignment",status:experimental,kind:static_,parameters:[{name:atom_indices,type:Option<Vec<usize>>,default:none},{name:conformer_ids,type:Option<Vec<usize>>,default:none},{name:weights,type:Option<Vec<f64>>,default:none},{name:reflect,type:bool,default:boolean(false)},{name:max_iterations,type:u32,default:integer(50)}],output:crate::ConformerAlignmentParameters,error:none,state:value_returning,operation:none,signature:fn(Option<Vec<usize>>,Option<Vec<usize>>,Option<Vec<f64>>,bool,u32)->crate::ConformerAlignmentParameters,},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.MaccsFingerprintParams",item:type,owner:type_,rust:crate::MaccsFingerprintParams,python:"MaccsFingerprintParams",javascript:"MaccsFingerprintParams",feature:"cap-fingerprints",status:experimental,role:parameter,},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"types.MaccsFingerprintError",item:type,owner:type_,rust:crate::MaccsFingerprintError,python:"MaccsFingerprintError",javascript:"MaccsFingerprintError",feature:"cap-fingerprints",status:experimental,role:error,},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"Molecule.maccs_fingerprint",item:callable,owner:molecule,rust:crate::Molecule::maccs_fingerprint,python:"maccs_fingerprint",javascript:"maccsFingerprint",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:crate::Fingerprint,error:crate::MaccsFingerprintError,state:read_only,operation:none,signature:for<'a> fn(&'a crate::Molecule)->Result<crate::Fingerprint,crate::MaccsFingerprintError>,},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"Molecule.maccs_fingerprint_raw",item:callable,owner:molecule,rust:crate::Molecule::maccs_fingerprint_raw,python:"maccs_fingerprint_raw",javascript:"maccsFingerprintRaw",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:crate::Fingerprint,error:crate::MaccsFingerprintError,state:read_only,operation:none,signature:for<'a> fn(&'a crate::Molecule)->Result<crate::Fingerprint,crate::MaccsFingerprintError>,},
        #[cfg(feature = "cap-fingerprints")]
        {semantic_id:"Molecule.maccs_fingerprint_with_params",item:callable,owner:molecule,rust:crate::Molecule::maccs_fingerprint_with_params,python:"maccs_fingerprint_with_params",javascript:"maccsFingerprintWithParams",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:params,type:&crate::MaccsFingerprintParams,default:required}],output:crate::Fingerprint,error:crate::MaccsFingerprintError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::Molecule,&'b crate::MaccsFingerprintParams)->Result<crate::Fingerprint,crate::MaccsFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.LayeredFingerprintParams",item:type,owner:type_,rust:crate::LayeredFingerprintParams,python:"LayeredFingerprintParams",javascript:"LayeredFingerprintParams",feature:"cap-fingerprints",status:experimental,role:parameter,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.LayeredFingerprintLayers",item:type,owner:type_,rust:crate::LayeredFingerprintLayers,python:"LayeredFingerprintLayers",javascript:"LayeredFingerprintLayers",feature:"cap-fingerprints",status:experimental,role:value,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.LayeredFingerprintResult",item:type,owner:type_,rust:crate::LayeredFingerprintResult,python:"LayeredFingerprintResult",javascript:"LayeredFingerprintResult",feature:"cap-fingerprints",status:experimental,role:result,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"types.LayeredFingerprintError",item:type,owner:type_,rust:crate::LayeredFingerprintError,python:"LayeredFingerprintError",javascript:"LayeredFingerprintError",feature:"cap-fingerprints",status:experimental,role:error,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.layered_fingerprint",item:callable,owner:molecule,rust:crate::Molecule::layered_fingerprint,python:"layered_fingerprint",javascript:"layeredFingerprint",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:crate::Fingerprint,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a> fn(&'a crate::Molecule)->Result<crate::Fingerprint,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.layered_fingerprint_with_params",item:callable,owner:molecule,rust:crate::Molecule::layered_fingerprint_with_params,python:"layered_fingerprint_with_params",javascript:"layeredFingerprintWithParams",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:params,type:&crate::LayeredFingerprintParams,default:required}],output:crate::Fingerprint,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::Molecule,&'b crate::LayeredFingerprintParams)->Result<crate::Fingerprint,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.layered_fingerprint_with_output",item:callable,owner:molecule,rust:crate::Molecule::layered_fingerprint_with_output,python:"layered_fingerprint_with_output",javascript:"layeredFingerprintWithOutput",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:crate::LayeredFingerprintResult,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a> fn(&'a crate::Molecule)->Result<crate::LayeredFingerprintResult,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"Molecule.layered_fingerprint_with_output_with_params",item:callable,owner:molecule,rust:crate::Molecule::layered_fingerprint_with_output_with_params,python:"layered_fingerprint_with_output_with_params",javascript:"layeredFingerprintWithOutputWithParams",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[{name:params,type:&crate::LayeredFingerprintParams,default:required}],output:crate::LayeredFingerprintResult,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::Molecule,&'b crate::LayeredFingerprintParams)->Result<crate::LayeredFingerprintResult,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"layered_query_fingerprint_with_params",item:callable,owner:module,rust:crate::layered_query_fingerprint_with_params,python:"layered_query_fingerprint_with_params",javascript:"layeredQueryFingerprintWithParams",feature:"cap-fingerprints",status:experimental,kind:module,parameters:[{name:query,type:&crate::QueryGraph,default:required},{name:params,type:&crate::LayeredFingerprintParams,default:required}],output:crate::Fingerprint,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::QueryGraph,&'b crate::LayeredFingerprintParams)->Result<crate::Fingerprint,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"layered_query_fingerprint_with_output_with_params",item:callable,owner:module,rust:crate::layered_query_fingerprint_with_output_with_params,python:"layered_query_fingerprint_with_output_with_params",javascript:"layeredQueryFingerprintWithOutputWithParams",feature:"cap-fingerprints",status:experimental,kind:module,parameters:[{name:query,type:&crate::QueryGraph,default:required},{name:params,type:&crate::LayeredFingerprintParams,default:required}],output:crate::LayeredFingerprintResult,error:crate::LayeredFingerprintError,state:read_only,operation:none,signature:for<'a,'b> fn(&'a crate::QueryGraph,&'b crate::LayeredFingerprintParams)->Result<crate::LayeredFingerprintResult,crate::LayeredFingerprintError>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"LayeredFingerprintResult.fingerprint",item:callable,owner:type_,rust:crate::LayeredFingerprintResult::fingerprint,python:"fingerprint",javascript:"fingerprint",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:&'a crate::Fingerprint,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::LayeredFingerprintResult)->&'a crate::Fingerprint,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"LayeredFingerprintResult.atom_counts",item:callable,owner:type_,rust:crate::LayeredFingerprintResult::atom_counts,python:"atom_counts",javascript:"atomCounts",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:shared,parameters:[],output:Option<&'a [u32]>,error:none,state:read_only,operation:none,signature:for<'a> fn(&'a crate::LayeredFingerprintResult)->Option<&'a [u32]>,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"LayeredFingerprintLayers.bits",item:callable,owner:type_,rust:crate::LayeredFingerprintLayers::bits,python:"bits",javascript:"bits",feature:"cap-fingerprints",status:experimental,kind:instance,receiver:owned,parameters:[],output:u32,error:none,state:value_returning,operation:none,signature:fn(crate::LayeredFingerprintLayers)->u32,},
        #[cfg(feature="cap-fingerprints")]
        {semantic_id:"LayeredFingerprintLayers.from_bits_retain",item:callable,owner:type_,rust:crate::LayeredFingerprintLayers::from_bits_retain,python:"from_bits_retain",javascript:"fromBitsRetain",feature:"cap-fingerprints",status:experimental,kind:static_,parameters:[{name:bits,type:u32,default:required}],output:crate::LayeredFingerprintLayers,error:none,state:value_returning,operation:none,signature:fn(u32)->crate::LayeredFingerprintLayers,},
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Fingerprint.from_on_bits", item: callable, owner: type_,
            rust: crate::Fingerprint::from_on_bits, python: "from_on_bits", javascript: "fromOnBits",
            feature: "cap-fingerprints", status: experimental, kind: static_,
            parameters: [{ name: n_bits, type: u32, default: required }, { name: on_bits, type: Vec<u32>, default: required }],
            output: crate::Fingerprint, error: crate::FingerprintError,
            state: value_returning, operation: none,
            signature: fn(u32, Vec<u32>) -> Result<crate::Fingerprint, crate::FingerprintError>,
        },
        #[cfg(feature = "cap-fingerprints")]
        {
            semantic_id: "Fingerprint.tanimoto", item: callable, owner: type_,
            rust: crate::Fingerprint::tanimoto, python: "tanimoto", javascript: "tanimoto",
            feature: "cap-fingerprints", status: experimental, kind: instance, receiver: shared,
            parameters: [{ name: other, type: &crate::Fingerprint, default: required }],
            output: f64, error: crate::FingerprintError,
            state: read_only, operation: none,
            signature: fn(&crate::Fingerprint, &crate::Fingerprint) -> Result<f64, crate::FingerprintError>,
        },
    ];
}
