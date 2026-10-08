import {Molecule,MorganParams,MorganFingerprintParams,MorganCallParams,MorganFingerprintGenerator,MorganSettings,MorganAtomInvariantsGenerator,MorganBondInvariantsGenerator,FingerprintAdditionalOutput,Fingerprint,SparseBitFingerprint,SparseCountFingerprint,SparseCountFingerprint32} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new MorganFingerprintParams(),g=MorganFingerprintGenerator.new(),call=MorganCallParams.new(),out=FingerprintAdditionalOutput.new();new MorganFingerprintGenerator(new MorganParams(),MorganAtomInvariantsGenerator.connectivity(),MorganBondInvariantsGenerator.new());MorganFingerprintGenerator.fromJson(g.toJson());const info:string=g.infoString();const settings:MorganSettings=g.settings();
const morgan_fingerprint_with_generator:Fingerprint=m.morganFingerprintWithGenerator(g,call,out);
const morgan_count_fingerprint_with_generator:SparseCountFingerprint32=m.morganCountFingerprintWithGenerator(g,call,out);
const morgan_sparse_fingerprint_with_generator:SparseBitFingerprint=m.morganSparseFingerprintWithGenerator(g,call,out);
const morgan_sparse_count_fingerprint_with_generator:SparseCountFingerprint=m.morganSparseCountFingerprintWithGenerator(g,call,out);
const morgan_sparse_count_fingerprint:SparseCountFingerprint=m.morganSparseCountFingerprint();
const morgan_sparse_count_fingerprint_with_params:SparseCountFingerprint=m.morganSparseCountFingerprintWithParams(p,out);
const morgan_sparse_fingerprint:SparseBitFingerprint=m.morganSparseFingerprint();
const morgan_sparse_fingerprint_with_params:SparseBitFingerprint=m.morganSparseFingerprintWithParams(p,out);
const morgan_count_fingerprint:SparseCountFingerprint32=m.morganCountFingerprint();
const morgan_count_fingerprint_with_params:SparseCountFingerprint32=m.morganCountFingerprintWithParams(p,out);
const morgan_fingerprint:Fingerprint=m.morganFingerprint();
const morgan_fingerprint_with_params:Fingerprint=m.morganFingerprintWithParams(p,out);
const fingerprints:(Fingerprint|null)[]=g.fingerprints([m,null],1);
const sparseFingerprints:(SparseBitFingerprint|null)[]=g.sparseFingerprints([m,null],1);
const counts:(SparseCountFingerprint32|null)[]=g.counts([m,null],1);
const sparseCounts:(SparseCountFingerprint|null)[]=g.sparseCounts([m,null],1);
const radius:number=settings.radius;settings.setRadius(1);settings.radius=1;
const onlyNonzeroInvariants:boolean=settings.onlyNonzeroInvariants;settings.setOnlyNonzeroInvariants(true);settings.onlyNonzeroInvariants=true;
const includeRedundantEnvironments:boolean=settings.includeRedundantEnvironments;settings.setIncludeRedundantEnvironments(true);settings.includeRedundantEnvironments=true;
const includeChirality:boolean=settings.includeChirality;settings.setIncludeChirality(true);settings.includeChirality=true;
const countSimulation:boolean=settings.countSimulation;settings.setCountSimulation(true);settings.countSimulation=true;
const fpSize:number=settings.fpSize;settings.setFpSize(64);settings.fpSize=64;
const bitsPerFeature:number=settings.bitsPerFeature;settings.setBitsPerFeature(2);settings.bitsPerFeature=2;
const countBounds:number[]=settings.countBounds;settings.setCountBounds([1,2]);settings.countBounds=[1,2];
const snap:MorganParams=settings.params();
// @ts-expect-error scalar params output requires explicit null or output
m.morganFingerprintWithParams(p);
// @ts-expect-error typed nullable molecule rows
g.fingerprints([{}]);
// @ts-expect-error bool input rejects coercion
settings.setCountSimulation(1);
