import {BioStructure,BioPdbWriteParams,BioMmcifWriteParams,BioPdbWriteError,BioMmcifWriteError} from 'cosmolkit-generated';
const s=BioStructure.fromPdb('END'),p=new BioPdbWriteParams(false,false,true,true,false),c=new BioMmcifWriteParams(false,true,null);p.preserveSerial=true;c.alignPairs=65535;c.authAll=true;const pdb:string=s.toPdb(),explicit:string=s.toPdbWithParams(p),cif:string=s.toMmcif(),selected:string=s.toMmcifWithParams(c);const a:void=s.writePdb('a'),d:void=s.writePdbWithParams('b',p),e:void=s.writeMmcif('c'),f:void=s.writeMmcifWithParams('d',c);declare const pe:BioPdbWriteError;declare const ce:BioMmcifWriteError;const serial:number|null=pe.serial,field:string|null=pe.detail,path:string|null=ce.path,cause:Error|null=ce.cause;
// @ts-expect-error boolean setter
p.endRecord=1;
// @ts-expect-error bounded integer setter
c.alignLoops='1';
// @ts-expect-error nullable context
const bad:number=pe.serial;
void [pdb,explicit,cif,selected,a,d,e,f,serial,field,path,cause,bad];
