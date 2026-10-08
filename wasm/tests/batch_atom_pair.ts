import {MoleculeBatch,AtomPairParams,AtomPairFingerprintParams,AtomPairAtomInvariantsGenerator,BatchQueryParams,Fingerprint,SparseBitFingerprint,SparseCountFingerprint,SparseCountFingerprint32,BatchFingerprintOutput,BatchFingerprintAdditionalOutput,BatchFingerprintOutputError} from 'cosmolkit-generated';
const generator=new AtomPairParams(),inv=new AtomPairAtomInvariantsGenerator(true,false),options=new AtomPairFingerprintParams(generator,[0],[],null,[], -1,inv,true),query=new BatchQueryParams(),value=MoleculeBatch.fromSmilesList(['CCO']);
const dense:(Fingerprint|null)[][]=[value.fingerprintAtomPairList(),value.fingerprintAtomPairListWithParams(options,query)];
const counts:(SparseCountFingerprint|null)[][]=[value.fingerprintAtomPairSparseCountList(),value.fingerprintAtomPairSparseCountListWithParams(options,query)];
const shortCounts:(SparseCountFingerprint32|null)[][]=[value.fingerprintAtomPairCountList(),value.fingerprintAtomPairCountListWithParams(options,query)];
const sparse:(SparseBitFingerprint|null)[][]=[value.fingerprintAtomPairSparseBitsList(),value.fingerprintAtomPairSparseBitsListWithParams(options,query)];
const results:(BatchFingerprintOutput|null)[][]=[value.fingerprintAtomPairWithOutputList(),value.fingerprintAtomPairWithOutputListWithParams(options,true,query)];
declare const result:BatchFingerprintOutput;const fp:Fingerprint=result.fingerprint(),extra:BatchFingerprintAdditionalOutput=result.additionalOutput(),atoms:number[]|null=extra.atomCounts(),bits:bigint[][]|null=extra.atomToBits(),info:Map<bigint,[number,number][]>|null=extra.bitInfoMap(),paths:Map<bigint,number[][]>|null=extra.bitPaths(),per:Map<bigint,number[][]>|null=extra.atomsPerBit();
declare const error:BatchFingerprintOutputError;const kind:string=error.kind;
const clone:AtomPairParams=generator.withJson(generator.toJson()),from:number[]|null=options.fromAtoms;
// @ts-expect-error frozen parameters
generator.fpSize=3;
// @ts-expect-error boolean collection flag
value.fingerprintAtomPairWithOutputListWithParams(options,'false',query);
void [dense,counts,shortCounts,sparse,results,fp,atoms,bits,info,paths,per,kind,clone,from];
