import {MoleculeBatch,SmilesWriteParams,BatchQueryParams,DrawingError} from "cosmolkit-generated";
const write=new SmilesWriteParams(),query=new BatchQueryParams(),configured=new BatchQueryParams(2,true,()=>{});
const root:number|null=write.rootedAtAtom,callback:(()=>void)|null=configured.progressCallback,jobs:number|null=query.nJobs,progress:boolean|null=query.progressBar;
const value=MoleculeBatch.fromSmilesList(["C"]);
const texts:(string|null)[][]=[value.toSmilesList(),value.toSmilesListWithParams(write,query),value.toSvgList(120,100),value.toSvgListWithParams(120,100,query)];
const matrices:(number[][]|null)[][]=[value.dgBoundsMatrixList(),value.dgBoundsMatrixListWithParams(query)];
declare const error:DrawingError;const domain:string=error.domain,kind:string=error.kind,message:string=error.message,cause:Error|null=error.cause;
// @ts-expect-error frozen writer value
write.canonical=false;
// @ts-expect-error callbacks must be callable
new BatchQueryParams(null,null,{});
// @ts-expect-error dimensions must be numbers
value.toSvgList("120",100);
void [root,callback,jobs,progress,texts,matrices,domain,kind,message,cause];
