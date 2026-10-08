import {Molecule,Coordinate2DParams,Coordinate2DError,Coordinate2DLayoutError,Coordinate2DTemplateError,DrawingError,DrawingWriteError} from 'cosmolkit-generated';
const p=new Coordinate2DParams(new Map([[0,new Float64Array([1,2])]]),true,false,3,4,-5,true,true,true);const a:Map<number,number[]>=p.coordinateMap;const bool:boolean=p.forceRdkit;const n:number=p.sampleSeed;
const m=Molecule.fromSmiles('CCO');const x:Molecule=m.with2dCoordinates();const y:Molecule=m.with2dCoordinatesWithParams(p);const u:void=m.compute2dCoordinates();const v:void=m.compute2dCoordinatesWithParams(p);const svg:string=m.toSvg(120,100);const png:Uint8Array=m.toPng(120,100);m.writeSvg('a.svg',120,100);m.writePng('a.png',120,100);
declare const e:Coordinate2DLayoutError;const atom:number|null=e.atom;const count:number|null=e.atomCount;const key:string|null=e.key;const cause:Error|null=e.cause;declare const t:Coordinate2DTemplateError;const path:string|null=t.path;const row:string|null=t.row;const index:number|null=t.index;const line:number|null=t.line;
// @ts-expect-error frozen source parameters
p.samples=4;
// @ts-expect-error dimensions are required
m.toSvg();
// @ts-expect-error numeric width
m.toPng('120',100);
void [a,bool,n,x,y,u,v,svg,png,atom,count,key,cause,path,row,index,line,Coordinate2DError,DrawingError,DrawingWriteError];
