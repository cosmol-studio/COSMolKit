import {Molecule,Fingerprint,MaccsFingerprintParams,MaccsFingerprintError} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new MaccsFingerprintParams(),width:number=p.nBits;
const a:Fingerprint=m.fingerprintMaccs(),b:Fingerprint=m.fingerprintMaccsRaw(),c:Fingerprint=m.fingerprintMaccsWithParams(p);
declare const e:MaccsFingerprintError;const option:string|null=e.option,reason:string|null=e.reason,bit:number|null=e.bit,cause:Error|null=e.cause;
// @ts-expect-error immutable configuration
p.nBits=167;
// @ts-expect-error no string coercion
new MaccsFingerprintParams('166');
