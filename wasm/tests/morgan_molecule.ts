import {Molecule,MorganParams,MorganFingerprintParams,MorganCallParams,MorganFingerprintGenerator,MorganSettings,MorganAtomInvariantsGenerator,MorganBondInvariantsGenerator,FingerprintAdditionalOutput,Fingerprint,SparseBitFingerprint,SparseCountFingerprint,SparseCountFingerprint32} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new MorganFingerprintParams(),g=MorganFingerprintGenerator.new(),call=MorganCallParams.new(),out=FingerprintAdditionalOutput.new();new MorganFingerprintGenerator(new MorganParams(),MorganAtomInvariantsGenerator.connectivity(),MorganBondInvariantsGenerator.new());MorganFingerprintGenerator.fromJson(g.toJson());const info:string=g.infoString();const settings:MorganSettings=g.settings();
const fingerprint_morgan_with_generator:Fingerprint=m.fingerprintMorganWithGenerator(g,call,out);
const fingerprint_morgan_count_with_generator:SparseCountFingerprint32=m.fingerprintMorganCountWithGenerator(g,call,out);
const fingerprint_morgan_sparse_with_generator:SparseBitFingerprint=m.fingerprintMorganSparseWithGenerator(g,call,out);
const fingerprint_morgan_sparse_count_with_generator:SparseCountFingerprint=m.fingerprintMorganSparseCountWithGenerator(g,call,out);
const fingerprint_morgan_sparse_count:SparseCountFingerprint=m.fingerprintMorganSparseCount();
const fingerprint_morgan_sparse_count_with_params:SparseCountFingerprint=m.fingerprintMorganSparseCountWithParams(p,out);
const fingerprint_morgan_sparse:SparseBitFingerprint=m.fingerprintMorganSparse();
const fingerprint_morgan_sparse_with_params:SparseBitFingerprint=m.fingerprintMorganSparseWithParams(p,out);
const fingerprint_morgan_count:SparseCountFingerprint32=m.fingerprintMorganCount();
const fingerprint_morgan_count_with_params:SparseCountFingerprint32=m.fingerprintMorganCountWithParams(p,out);
const fingerprint_morgan:Fingerprint=m.fingerprintMorgan();
const fingerprint_morgan_with_params:Fingerprint=m.fingerprintMorganWithParams(p,out);
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
m.fingerprintMorganWithParams(p);
// @ts-expect-error typed nullable molecule rows
g.fingerprints([{}]);
// @ts-expect-error bool input rejects coercion
settings.setCountSimulation(1);
