import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('All dataset/stream opens preserve real WASM file errors and reusable read policies',()=>{
 const path='/supplier-input.sdf',p=new b.SdfReadParams(false,false,false,true,false,b.SdfCoordinateMode.Require3D);
 for(const cls of [b.SdfDataset,b.SdfRecordStream])for(const fn of [()=>cls.open(path),()=>cls.openWithParams(path,p),()=>cls.openWithParams(path,p)]){
 let e;try{fn();}catch(x){e=x;}assert.ok(e instanceof Error);assert.equal(e.name,'MolecularIoError');assert.equal(e.domain,'io');assert.equal(e.kind,'Io');assert.ok(e.detail instanceof b.MolecularIoError);assert.equal(e.filename,path);assert.equal(e.detail.filename,path);assert.equal(e.detail.ioKind,'Unsupported');assert.equal(e.detail.errno,null);assert.ok(e.cause instanceof Error);
 }
 assert.equal(p.sanitize,false);assert.equal(p.coordinateMode,b.SdfCoordinateMode.Require3D);
});
test('Reusable reader stores source path and full detached params without touching filesystem',()=>{
 const path='/missing-is-allowed.sdf',defaults=b.SdfReader.open(path);assert.equal(defaults.path(),path);assert.equal(defaults.params().sanitize,true);assert.equal(defaults.params().coordinateMode,b.SdfCoordinateMode.Preserve);
 const p=new b.SdfReadParams(false,false,false,true,false,b.SdfCoordinateMode.Require3D),reader=b.SdfReader.openWithParams(path,p);assert.equal(reader.path(),path);
 for(let repeat=0;repeat<2;repeat++){const value=reader.params();assert.equal(value.sanitize,false);assert.equal(value.removeHydrogens,false);assert.equal(value.strictParsing,false);assert.equal(value.expandAttachmentPoints,true);assert.equal(value.processPropertyLists,false);assert.equal(value.coordinateMode,b.SdfCoordinateMode.Require3D);value.free();}
 p.free();assert.equal(reader.params().coordinateMode,b.SdfCoordinateMode.Require3D);assert.equal(reader.path(),path);
});
test('Supplier exports retain complete registry methods and class boundary validation',()=>{
 const methods={SdfDataset:['len','isEmpty','path','iter','metadata','record','recordWithParams','recordText'],SdfRecordStream:['nextRecord','isEnd','recordsConsumed','bytesConsumed','linesConsumed'],SdfReader:['path','params'],SdfRecordMetadata:['index','byteOffset','byteLen','byteRange','lineRange','title'],SdfDatasetIterator:['next']};
 for(const [cls,names]of Object.entries(methods)){assert.equal(typeof b[cls],'function');for(const name of names)assert.equal(typeof b[cls].prototype[name],'function',`${cls}.${name}`);}
 for(const cls of [b.SdfDataset,b.SdfRecordStream,b.SdfReader])for(const v of [null,undefined,{},new b.XyzWriteParams()])assert.throws(()=>cls.openWithParams('/input.sdf',v),TypeError);
});
