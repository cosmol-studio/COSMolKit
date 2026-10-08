import {MoleculeBatch,BatchParams,SanitizeParams,SanitizeOperations,SanitizeStage,AddHsParams,RemoveHsParams,KekulizeParams,Coordinate2DParams,SanitizeError,HydrogenError,KekulizeError,Coordinate2DError} from "cosmolkit-generated";
const operations:SanitizeOperations=SanitizeOperations.CLEANUP.or(SanitizeOperations.PROPERTIES);
const bits:number=operations.bits();const empty:boolean=operations.isEmpty();const subset:boolean=operations.contains(SanitizeOperations.ALL);
const sanitize=new SanitizeParams(operations), add=new AddHsParams(false,false,false,false,new Uint32Array([0])),remove=new RemoveHsParams(),kek=new KekulizeParams(true,true,100),coords=new Coordinate2DParams(new Map([[0,[1,2]]]),false,true,0,0,-1);
const atomIds:number[]|null=add.onlyOnAtoms;const map:Map<number,[number,number]>=coords.coordinateMap;const stage:SanitizeStage=SanitizeStage.Properties;
const batch=MoleculeBatch.fromSmilesList(["C"]),params=new BatchParams();
const results:MoleculeBatch[]=[batch.sanitize(),batch.sanitizeWithParams(sanitize,params),batch.withHydrogens(),batch.withHydrogensWithParams(add,params),batch.withoutHydrogens(),batch.withoutHydrogensWithParams(remove,params),batch.withKekulizedBonds(),batch.withKekulizedBondsWithParams(kek,params),batch.with2dCoordinates(),batch.with2dCoordinatesWithParams(coords,params)];
declare const errors: [SanitizeError,HydrogenError,KekulizeError,Coordinate2DError];for(const error of errors){const domain:string=error.domain,kind:string=error.kind,message:string=error.message,cause:Error|null=error.cause;void [domain,kind,message,cause];}
// @ts-expect-error frozen parameter
coords.samples=2;
// @ts-expect-error frozen parameter
remove.sanitize=false;
// @ts-expect-error integer input must be a number
new KekulizeParams(true,true,"100");
// @ts-expect-error coordinate map requires two-dimensional tuples
new Coordinate2DParams(new Map([[0,[1,2,3]]]));
void [bits,empty,subset,atomIds,map,stage,results];
