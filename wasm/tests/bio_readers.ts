import {BioStructure,BioReadParams,BioPdbReadParams,BioCoordinateFormat,BioReadError,BioPdbReadError,BioMmcifReadError,BioPdbReadStage,BioMmcifReadStage} from 'cosmolkit-generated';
const p=new BioPdbReadParams(0,false,false,false,false),r=new BioReadParams(BioCoordinateFormat.Mmcif,'named.cif');const values:BioStructure[]=[BioStructure.fromPdb('END'),BioStructure.fromPdbWithParams('END',p),BioStructure.fromMmcif('data_x'),BioStructure.fromText('END'),BioStructure.fromTextWithParams('data_x',r),BioStructure.read('a.pdb'),BioStructure.readWithFormat('a.pdb',BioCoordinateFormat.Pdb)];const f:BioCoordinateFormat=values[0].inputFormat(),name:string=values[0].name(),counts:number[]=[values[0].numAtoms(),values[0].numModels(),values[0].numChains(),values[0].numResidues(),values[0].numEntities()];declare const e:BioPdbReadError;const line:number|null=e.lineNumber(),tag:Uint8Array|null=e.recordTag(),stage:BioPdbReadStage=e.stage();declare const m:BioMmcifReadError;const other:BioMmcifReadStage=m.stage();declare const read:BioReadError;const kind:string=read.kind;
// @ts-expect-error readonly read options
r.sourceName='other';
// @ts-expect-error nullable record bytes
const bad:Uint8Array=e.recordTag();
// @ts-expect-error boolean flag type
new BioPdbReadParams(0,'false');
void [f,name,counts,line,tag,stage,other,kind,bad];
