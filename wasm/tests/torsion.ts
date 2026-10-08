import {Molecule,TopologicalTorsionParams,TopologicalTorsionFingerprintParams,TopologicalTorsionCallParams,TopologicalTorsionFingerprintGenerator,TopologicalTorsionSettings,AtomPairAtomInvariantsGenerator,LegacyTopologicalTorsionParams,SparseCountFingerprint,FingerprintAdditionalOutput,Fingerprint,SparseBitFingerprint,SparseCountFingerprint32} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new TopologicalTorsionFingerprintParams(),g=TopologicalTorsionFingerprintGenerator.new(),call=new TopologicalTorsionCallParams(),out=FingerprintAdditionalOutput.new();new TopologicalTorsionFingerprintGenerator(new TopologicalTorsionParams(),new AtomPairAtomInvariantsGenerator());TopologicalTorsionFingerprintGenerator.fromJson(g.toJson());const info:string=g.infoString();const settings:TopologicalTorsionSettings=g.settings();
const topological_torsion_fingerprint_with_generator:Fingerprint=m.topologicalTorsionFingerprintWithGenerator(g,call,out);
const topological_torsion_count_fingerprint_with_generator:SparseCountFingerprint32=m.topologicalTorsionCountFingerprintWithGenerator(g,call,out);
const topological_torsion_sparse_fingerprint_with_generator:SparseBitFingerprint=m.topologicalTorsionSparseFingerprintWithGenerator(g,call,out);
const topological_torsion_sparse_count_fingerprint_with_generator:SparseCountFingerprint=m.topologicalTorsionSparseCountFingerprintWithGenerator(g,call,out);
const topological_torsion_sparse_count_fingerprint:SparseCountFingerprint=m.topologicalTorsionSparseCountFingerprint();
const topological_torsion_sparse_count_fingerprint_with_params:SparseCountFingerprint=m.topologicalTorsionSparseCountFingerprintWithParams(p,out);
const topological_torsion_sparse_fingerprint:SparseBitFingerprint=m.topologicalTorsionSparseFingerprint();
const topological_torsion_sparse_fingerprint_with_params:SparseBitFingerprint=m.topologicalTorsionSparseFingerprintWithParams(p,out);
const topological_torsion_count_fingerprint:SparseCountFingerprint32=m.topologicalTorsionCountFingerprint();
const topological_torsion_count_fingerprint_with_params:SparseCountFingerprint32=m.topologicalTorsionCountFingerprintWithParams(p,out);
const topological_torsion_fingerprint:Fingerprint=m.topologicalTorsionFingerprint();
const topological_torsion_fingerprint_with_params:Fingerprint=m.topologicalTorsionFingerprintWithParams(p,out);
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
m.topologicalTorsionFingerprintWithParams(p);
// @ts-expect-error typed nullable molecule rows
g.fingerprints([{}]);
// @ts-expect-error bool input rejects coercion
settings.setCountSimulation(1);

const legacy=new LegacyTopologicalTorsionParams(),legacyStatic=LegacyTopologicalTorsionParams.new(3,true,128,2,[0],[],[1,2,3]);
const legacyDense:Fingerprint=m.legacyTopologicalTorsionFingerprint();const legacyDenseParams:Fingerprint=m.legacyTopologicalTorsionFingerprintWithParams(legacy);
const legacyCount:SparseCountFingerprint=m.legacyTopologicalTorsionCountFingerprint();const legacyCountParams:SparseCountFingerprint=m.legacyTopologicalTorsionCountFingerprintWithParams(legacy);
const legacySparse:SparseCountFingerprint=m.legacyTopologicalTorsionSparseCountFingerprint();const legacySparseParams:SparseCountFingerprint=m.legacyTopologicalTorsionSparseCountFingerprintWithParams(legacy);
const ids:bigint[]=m.topologicalTorsionIds();const sizedIds:bigint[]=m.topologicalTorsionIdsWithParams(3);
// @ts-expect-error path IDs use exact unsigned widths
m.topologicalTorsionIdsWithParams(3n);
