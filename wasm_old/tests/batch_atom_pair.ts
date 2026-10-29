import {MoleculeBatch,AtomPairGeneratorParams,AtomPairFingerprintParams,AtomPairAtomInvariantsGenerator,BatchQueryParams,Fingerprint,SparseBitFingerprint,SparseCountFingerprint,SparseCountFingerprint32,BatchFingerprintOutput,BatchFingerprintAdditionalOutput,BatchFingerprintOutputError} from 'cosmolkit-generated';
const generator=new AtomPairGeneratorParams(),inv=new AtomPairAtomInvariantsGenerator(true,false),options=new AtomPairFingerprintParams(generator,[0],[],null,[], -1,inv,true),query=new BatchQueryParams(),value=MoleculeBatch.fromSmilesList(['CCO']);
const dense:(Fingerprint|null)[][]=[value.fingerprintAtomPair(),value.fingerprintAtomPairWithParams(options,query)];
const counts:(SparseCountFingerprint|null)[][]=[value.fingerprintAtomPairSparseCount(),value.fingerprintAtomPairSparseCountWithParams(options,query)];
const shortCounts:(SparseCountFingerprint32|null)[][]=[value.fingerprintAtomPairCount(),value.fingerprintAtomPairCountWithParams(options,query)];
const sparse:(SparseBitFingerprint|null)[][]=[value.fingerprintAtomPairSparse(),value.fingerprintAtomPairSparseWithParams(options,query)];
const results:(BatchFingerprintOutput|null)[][]=[value.fingerprintAtomPairWithOutput(),value.fingerprintAtomPairWithOutputWithParams(options,true,query)];
declare const result:BatchFingerprintOutput;const fp:Fingerprint=result.fingerprint(),extra:BatchFingerprintAdditionalOutput=result.additionalOutput(),atoms:number[]|null=extra.atomCounts(),bits:bigint[][]|null=extra.atomToBits(),info:Map<bigint,[number,number][]>|null=extra.bitInfoMap(),paths:Map<bigint,number[][]>|null=extra.bitPaths(),per:Map<bigint,number[][]>|null=extra.atomsPerBit();
declare const error:BatchFingerprintOutputError;const kind:string=error.kind;
const clone:AtomPairGeneratorParams=generator.withJson(generator.toJson()),from:number[]|null=options.fromAtoms;
generator.fpSize=3;
// @ts-expect-error boolean collection flag
value.fingerprintAtomPairWithOutputWithParams(options,'false',query);
void [dense,counts,shortCounts,sparse,results,fp,atoms,bits,info,paths,per,kind,clone,from];
