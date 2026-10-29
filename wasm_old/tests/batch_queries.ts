import {MoleculeBatch,SmilesWriteParams,BatchQueryParams,DrawingError} from "cosmolkit-generated";
const write=new SmilesWriteParams(),query=new BatchQueryParams(),configured=new BatchQueryParams(2,true,()=>{});
const root:number|null=write.rootedAtAtom,callback:(()=>void)|null=configured.progressCallback,jobs:number|null=query.nJobs,progress:boolean|null=query.progressBar;
const value=MoleculeBatch.fromSmilesList(["C"]);
const texts:(string|null)[][]=[value.toSmiles(),value.toSmilesWithParams(write,query),value.toSvg(120,100),value.toSvgWithParams(120,100,query)];
const matrices:(number[][]|null)[][]=[value.dgBoundsMatrix(),value.dgBoundsMatrixWithParams(query)];
declare const error:DrawingError;const domain:string=error.domain,kind:string=error.kind,message:string=error.message,cause:Error|null=error.cause;
write.canonical=false;
// @ts-expect-error callbacks must be callable
new BatchQueryParams(null,null,{});
// @ts-expect-error dimensions must be numbers
value.toSvg("120",100);
void [root,callback,jobs,progress,texts,matrices,domain,kind,message,cause];
