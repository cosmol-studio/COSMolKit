import {Molecule,Mol2ReadParams,Mol2Type,XyzWriteParams,MolecularIoError,XyzReadError,XyzWriteError,Mol2ReadError,Mol2PostError} from 'cosmolkit-generated';
const r=new Mol2ReadParams(false,false,Mol2Type.Corina,false),w=new XyzWriteParams(null,2),id:number|null=w.conformerId,precision:number=w.precision,variant:Mol2Type=r.variant,sanitize:boolean=r.sanitize,remove:boolean=r.removeHs,cleanup:boolean=r.cleanupSubstructures;
const a:Molecule=Molecule.fromXyzBlock(''),b:Molecule=Molecule.readXyz(''),c:Molecule=Molecule.fromMol2(''),d:Molecule=Molecule.readMol2(''),e:Molecule=Molecule.fromMol2WithParams('',r),f:Molecule=Molecule.readMol2WithParams('',r);const text:string=a.toXyz(),custom:string=a.toXyzWithParams(w),write:void=a.writeXyz(''),writeCustom:void=a.writeXyzWithParams('',w);
declare const xr:XyzReadError,xw:XyzWriteError,mr:Mol2ReadError,mp:Mol2PostError,io:MolecularIoError;const line:number|null=xr.line,value:string|null=xr.value,missing:number|null=xw.id,feature:string|null=mr.feature,parse:string|null=mr.parseDetail,stage:string|null=mp.stage;for(const err of [xr,xw,mr,mp,io]){const domain:string=err.domain,kind:string=err.kind,message:string=err.message,cause:Error|null=err.cause;}
// @ts-expect-error typed class required
Molecule.fromMol2WithParams('',{});
// @ts-expect-error typed class required
a.toXyzWithParams({});
// @ts-expect-error usize is a checked number
new XyzWriteParams(0n);
// @ts-expect-error nullable selector
const selected:number=w.conformerId;
