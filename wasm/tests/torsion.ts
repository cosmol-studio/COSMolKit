import {Molecule,TopologicalTorsionParams,TopologicalTorsionFingerprintParams,TopologicalTorsionCallParams,TopologicalTorsionFingerprintGenerator,TopologicalTorsionSettings,AtomPairAtomInvariantsGenerator,LegacyTopologicalTorsionParams,SparseCountFingerprint,FingerprintAdditionalOutput,Fingerprint,SparseBitFingerprint,SparseCountFingerprint32} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new TopologicalTorsionFingerprintParams(),g=TopologicalTorsionFingerprintGenerator.new(),call=new TopologicalTorsionCallParams(),out=FingerprintAdditionalOutput.new();new TopologicalTorsionFingerprintGenerator(new TopologicalTorsionParams(),new AtomPairAtomInvariantsGenerator());TopologicalTorsionFingerprintGenerator.fromJson(g.toJson());const info:string=g.infoString();const settings:TopologicalTorsionSettings=g.settings();
const fingerprint_topological_torsion_with_generator:Fingerprint=m.fingerprintTopologicalTorsionWithGenerator(g,call,out);
const fingerprint_topological_torsion_count_with_generator:SparseCountFingerprint32=m.fingerprintTopologicalTorsionCountWithGenerator(g,call,out);
const fingerprint_topological_torsion_sparse_with_generator:SparseBitFingerprint=m.fingerprintTopologicalTorsionSparseWithGenerator(g,call,out);
const fingerprint_topological_torsion_sparse_count_with_generator:SparseCountFingerprint=m.fingerprintTopologicalTorsionSparseCountWithGenerator(g,call,out);
const fingerprint_topological_torsion_sparse_count:SparseCountFingerprint=m.fingerprintTopologicalTorsionSparseCount();
const fingerprint_topological_torsion_sparse_count_with_params:SparseCountFingerprint=m.fingerprintTopologicalTorsionSparseCountWithParams(p,out);
const fingerprint_topological_torsion_sparse:SparseBitFingerprint=m.fingerprintTopologicalTorsionSparse();
const fingerprint_topological_torsion_sparse_with_params:SparseBitFingerprint=m.fingerprintTopologicalTorsionSparseWithParams(p,out);
const fingerprint_topological_torsion_count:SparseCountFingerprint32=m.fingerprintTopologicalTorsionCount();
const fingerprint_topological_torsion_count_with_params:SparseCountFingerprint32=m.fingerprintTopologicalTorsionCountWithParams(p,out);
const fingerprint_topological_torsion:Fingerprint=m.fingerprintTopologicalTorsion();
const fingerprint_topological_torsion_with_params:Fingerprint=m.fingerprintTopologicalTorsionWithParams(p,out);
const fingerprints:(Fingerprint|null)[]=g.fingerprints([m,null],1);
const sparseFingerprints:(SparseBitFingerprint|null)[]=g.sparseFingerprints([m,null],1);
const counts:(SparseCountFingerprint32|null)[]=g.counts([m,null],1);
const sparseCounts:(SparseCountFingerprint|null)[]=g.sparseCounts([m,null],1);
const torsionAtomCount:number=settings.torsionAtomCount;settings.setTorsionAtomCount(3);settings.torsionAtomCount=4;
const onlyShortestPaths:boolean=settings.onlyShortestPaths;settings.setOnlyShortestPaths(true);settings.onlyShortestPaths=false;
const includeChirality:boolean=settings.includeChirality;settings.setIncludeChirality(true);settings.includeChirality=true;
const countSimulation:boolean=settings.countSimulation;settings.setCountSimulation(true);settings.countSimulation=true;
const fpSize:number=settings.fpSize;settings.setFpSize(64);settings.fpSize=64;
const bitsPerFeature:number=settings.bitsPerFeature;settings.setBitsPerFeature(2);settings.bitsPerFeature=2;
const countBounds:number[]=settings.countBounds;settings.setCountBounds([1,2]);settings.countBounds=[1,2];
const snap:TopologicalTorsionParams=settings.params();
// @ts-expect-error scalar params output requires explicit null or output
m.fingerprintTopologicalTorsionWithParams(p);
// @ts-expect-error typed nullable molecule rows
g.fingerprints([{}]);
// @ts-expect-error bool input rejects coercion
settings.setCountSimulation(1);

const legacy=new LegacyTopologicalTorsionParams(),legacyStatic=LegacyTopologicalTorsionParams.new(3,true,128,2,[0],[],[1,2,3]);
const legacyDense:Fingerprint=m.fingerprintTopologicalTorsionLegacy();const legacyDenseParams:Fingerprint=m.fingerprintTopologicalTorsionLegacyWithParams(legacy);
const legacyCount:SparseCountFingerprint=m.fingerprintTopologicalTorsionCountLegacy();const legacyCountParams:SparseCountFingerprint=m.fingerprintTopologicalTorsionCountLegacyWithParams(legacy);
const legacySparse:SparseCountFingerprint=m.fingerprintTopologicalTorsionSparseCountLegacy();const legacySparseParams:SparseCountFingerprint=m.fingerprintTopologicalTorsionSparseCountLegacyWithParams(legacy);
const ids:bigint[]=m.topologicalTorsionIds();const sizedIds:bigint[]=m.topologicalTorsionIdsWithParams(3);
// @ts-expect-error path IDs use exact unsigned widths
m.topologicalTorsionIdsWithParams(3n);
