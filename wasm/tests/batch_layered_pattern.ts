import {MoleculeBatch,LayeredFingerprintLayers,LayeredFingerprintParams,LayeredFingerprintResult,PatternFingerprintParams,LayeredFingerprintError,PatternFingerprintError,BatchQueryParams,Fingerprint} from 'cosmolkit-generated';
const mask=Fingerprint.fromOnBits(2048,[674]),options=new LayeredFingerprintParams(4294967295,1,7,2048,[0,0,0],mask,true,[]),query=new BatchQueryParams(),batch=MoleculeBatch.fromSmilesList(['CCO']),pattern=new PatternFingerprintParams(128,true);
const dense:(Fingerprint|null)[][]=[batch.fingerprintLayeredList(),batch.fingerprintLayeredListWithParams(options,query),batch.fingerprintPatternList(),batch.fingerprintPatternListWithParams(pattern,query)];
const outputs:(LayeredFingerprintResult|null)[][]=[batch.fingerprintLayeredWithOutputList(),batch.fingerprintLayeredWithOutputListWithParams(options,query)];
declare const result:LayeredFingerprintResult;const counts:number[]|null=result.atomCounts(),fp:Fingerprint=result.fingerprint(),flags:number=LayeredFingerprintLayers.fromBitsRetain(0xffffffff).bits(),borrowed:Fingerprint|null=options.setOnlyBits;
declare const e:LayeredFingerprintError;declare const pe:PatternFingerprintError;const kinds:string[]=[e.kind,pe.kind];
// @ts-expect-error readonly parameters
options.atomCounts=[];
// @ts-expect-error result nullable collection
const wrong:number[]=result.atomCounts();
// @ts-expect-error checked unsigned number
new PatternFingerprintParams(64n);
void [dense,outputs,counts,fp,flags,borrowed,kinds,wrong];
