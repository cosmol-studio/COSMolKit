import {MoleculeBatch,BatchImageParams,BatchExportReport,BatchParams,BatchError,BatchImageError,DrawingWriteError} from "cosmolkit-generated";
const defaults=new BatchImageParams(),options=new BatchImageParams("svg",120,100,new BatchParams(),["ethanol",null],"counts.json");
const names:(string|null)[]|null=options.filenames,path:string|null=options.reportPath,execution:BatchParams=options.execution;
const batch=MoleculeBatch.fromSmilesList(["C"]),report:BatchExportReport=batch.toImages("images"),explicit:BatchExportReport=batch.toImagesWithParams("images",options);
const counts:number[]=[report.total(),report.success(),report.failed(),report.written,report.skipped];const errors:BatchError[]=report.errors();const written:void=report.writeReport("counts.json");
declare const imageError:BatchImageError,writeError:DrawingWriteError;const causes:(Error|null)[]=[imageError.cause,writeError.cause];
// @ts-expect-error frozen values
options.width=7;
// @ts-expect-error filenames must be nullable strings
new BatchImageParams("png",300,300,null,[7]);
void [defaults,names,path,execution,explicit,counts,errors,written,causes];
