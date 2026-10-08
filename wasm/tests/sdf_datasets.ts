import {SdfDataset,SdfDatasetIterator,SdfRecordMetadata,SdfRecordStream,SdfReader,SdfRecord,SdfReadParams} from 'cosmolkit-generated';
const p=new SdfReadParams(),dataset:SdfDataset=SdfDataset.open(''),custom:SdfDataset=SdfDataset.openWithParams('',p),len:number=dataset.len(),empty:boolean=dataset.isEmpty(),path:string=dataset.path(),metadata:SdfRecordMetadata|null=dataset.metadata(0),record:SdfRecord=dataset.record(0),recordCustom:SdfRecord=dataset.recordWithParams(0,p),text:string=dataset.recordText(0),iterator:SdfDatasetIterator=dataset.iter(),next:IteratorResult<SdfRecord>=iterator.next();for(const record of iterator){const value:SdfRecord=record;}
declare const m:SdfRecordMetadata;const index:number=m.index(),offset:bigint=m.byteOffset(),length:bigint=m.byteLen(),bytes:[bigint,bigint]=m.byteRange(),lines:[number,number]=m.lineRange(),title:string|null=m.title();const stream:SdfRecordStream=SdfRecordStream.open(''),streamCustom:SdfRecordStream=SdfRecordStream.openWithParams('',p),item:SdfRecord|null=stream.nextRecord(),end:boolean=stream.isEnd(),records:number=stream.recordsConsumed(),bytesConsumed:bigint=stream.bytesConsumed(),linesConsumed:number=stream.linesConsumed();const reader:SdfReader=SdfReader.open(''),readerCustom:SdfReader=SdfReader.openWithParams('',p),readerPath:string=reader.path(),readerParams:SdfReadParams=reader.params();
// @ts-expect-error u64 metadata stays bigint
const lossy:number=m.byteOffset();
// @ts-expect-error dataset metadata can be absent
const mandatory:SdfRecordMetadata=dataset.metadata(0);
// @ts-expect-error stream completion is nullable
const required:SdfRecord=stream.nextRecord();
// @ts-expect-error typed read policy required
SdfDataset.openWithParams('',{});
// @ts-expect-error usize index is a number
dataset.record(0n);
