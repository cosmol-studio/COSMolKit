import{MoleculeBatch,BatchRecord,BatchErrorMode,BatchExportParams,BatchExportReport,SdfFormat,SdfReadParams,SdfDataset,SdfReader,SdfRecordStream,SdfBatchIterator,SdfReaderBatchIterator}from'cosmolkit-generated';
const p=new SdfReadParams(),params=new BatchExportParams(SdfFormat.V3000,BatchErrorMode.KeepErrors,null,false),format:SdfFormat=params.format,mode:BatchErrorMode=params.errors,n:number|null=params.nJobs,progress:boolean|null=params.progressBar,a:MoleculeBatch=MoleculeBatch.fromSdfRecords(''),b:MoleculeBatch=MoleculeBatch.fromSdfRecordsWithParams('',p,BatchErrorMode.Strict,null),c:MoleculeBatch=MoleculeBatch.readSdf(''),d:MoleculeBatch=MoleculeBatch.readSdfWithParams('',p,BatchErrorMode.Strict,null,false);declare const ds:SdfDataset;const selected:MoleculeBatch=MoleculeBatch.fromDatasetIndices(ds,new Uint32Array([0]),BatchErrorMode.Strict),record:BatchRecord|null=a.get(0),records:BatchRecord[]=a.records();const r:BatchExportReport=a.toSdf(''),rf:BatchExportReport=a.toSdfFiles(''),rp:BatchExportReport=a.toSdfWithParams('',params,null),rfp:BatchExportReport=a.toSdfFilesWithParams('',params,['x',null],null),total:number=r.total(),success:number=r.success(),failed:number=r.failed();const it:SdfBatchIterator=ds.batches(2,BatchErrorMode.Strict,null,false),item:MoleculeBatch|null=it.nextBatch(),result:IteratorResult<MoleculeBatch>=it.next();for(const item of it){const v:MoleculeBatch=item;}const reader=SdfReader.open(''),defaultIterator:SdfReaderBatchIterator=reader.batches(),customIterator:SdfReaderBatchIterator=reader.batches(2,BatchErrorMode.KeepErrors,null);declare const stream:SdfRecordStream;const streamIterator:SdfReaderBatchIterator=stream.batches(2,BatchErrorMode.Strict,null),next:MoleculeBatch|null=streamIterator.nextBatch();
// @ts-expect-error nullable batch record
const certain:BatchRecord=a.get(0);
// @ts-expect-error typed export policy required
a.toSdfWithParams('',{},null);
// @ts-expect-error required nullable report argument
a.toSdfWithParams('',params);
// @ts-expect-error typed reader policy required
MoleculeBatch.fromSdfRecordsWithParams('',{},BatchErrorMode.Strict,null);
// @ts-expect-error unsigned indices use numbers
a.get(0n);
