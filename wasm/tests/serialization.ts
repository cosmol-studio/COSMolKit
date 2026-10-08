import {Molecule,PickleError} from 'cosmolkit-generated';
const m:Molecule=Molecule.fromSmiles('CCO'),data:Uint8Array=m.toBinary(),copy:Molecule=Molecule.fromBinary(data),slice:Molecule=Molecule.fromBinary(data.subarray(0));
declare const e:PickleError;
const domain:string=e.domain,kind:string=e.kind,message:string=e.message,cause:Error|null=e.cause,version:number|null=e.version,major:number|null=e.major,minor:number|null=e.minor,section:number|null=e.section,expected:number|null=e.expected,actual:number|null=e.actual,value:number|null=e.value,typeName:string|null=e.typeName,count:number|null=e.count;
// @ts-expect-error Uint8Array required
Molecule.fromBinary([1,2]);
// @ts-expect-error bytes required
Molecule.fromBinary();
// @ts-expect-error immutable error payload
 e.version=1;
// @ts-expect-error return is a byte array
const text:string=m.toBinary();
