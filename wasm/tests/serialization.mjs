import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
function failure(data,kind){let e;try{b.Molecule.fromBinary(data);}catch(x){e=x;}assert.ok(e instanceof Error);assert.equal(e.name,'PickleError');assert.equal(e.domain,'serialization');if(kind)assert.equal(e.kind,kind);assert.ok(e.detail instanceof b.PickleError);assert.equal(e.detail.kind,e.kind);assert.equal(e.detail.message,e.message);assert.equal(e.detail.cause,null);return e;}
test('Binary round trips preserve stereo, isotope, mapping, ring caches, independent bytes and object lifetime',()=>{
 for(const smiles of ['', 'CCO','c1ccccc1','F[C@H](Cl)Br','[13CH3:7][NH3+]','C.O']){
  const m=b.Molecule.fromSmiles(smiles),before=m.toSmiles(),bytes=m.toBinary();assert.ok(bytes instanceof Uint8Array);assert.deepEqual([...bytes.slice(0,12)],[67,79,83,77,79,76,0,0,2,0,0,0]);
  const restored=b.Molecule.fromBinary(bytes);assert.equal(restored.toSmiles(),before);assert.equal(restored.numAtoms(),m.numAtoms());assert.equal(restored.numBonds(),m.numBonds());assert.equal(restored.numRings(),m.numRings());assert.deepEqual(restored.toBinary(),bytes);
  const another=m.toBinary();another[0]=255;assert.deepEqual(m.toBinary(),bytes);assert.equal(m.toSmiles(),before);m.free();assert.equal(restored.toSmiles(),before);const again=b.Molecule.fromBinary(bytes);bytes.fill(0);assert.equal(again.toSmiles(),before);
 }
});
test('Binary preserves exact 3D coordinates, ordered duplicate SDF fields and serialized bytes',()=>{
 const m=b.Molecule.fromXyzBlock('2\narchive\nC -0.0 1.25 -2.5\nO 3 4 5\n'),before=Array.from(m.coordinates3d(0)),bytes=m.toBinary(),restored=b.Molecule.fromBinary(bytes);assert.deepEqual(Array.from(restored.coordinates3d(0)),before);assert.deepEqual(restored.toBinary(),bytes);assert.deepEqual(Array.from(m.coordinates3d(0)),before);
 const sdf='binary\n  COSMolKit         2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <NOTE>\nfirst\n\n>  <NOTE>\nsecond\n\n$$$$\n';const source=b.Molecule.fromSdf(sdf),binary=source.toBinary(),copy=b.Molecule.fromBinary(binary);assert.equal(copy.toSdf(),source.toSdf());assert.deepEqual(b.SdfRecord.fromSdf(copy.toSdf()).dataFields(),[['NOTE','first'],['NOTE','second']]);assert.deepEqual(copy.toBinary(),binary);
});
test('Malformed archives preserve exact error classification and payloads without changing existing values',()=>{
 const m=b.Molecule.fromSmiles('CCO'),bytes=m.toBinary();failure(new Uint8Array(),'UnexpectedEof');const legacy=failure(Uint8Array.of(255),'UnsupportedVersion');assert.equal(legacy.detail.version,255);assert.equal(legacy.version,255);assert.equal(legacy.detail.major,null);assert.equal(legacy.detail.typeName,null);assert.equal(legacy.detail.count,null);
 for(let end=0;end<bytes.length;end++)failure(bytes.slice(0,end));
 let data=bytes.slice();data[8]=99;const archive=failure(data,'UnsupportedArchiveVersion');assert.equal(archive.detail.major,99);assert.equal(archive.detail.minor,0);
 data=bytes.slice();data[16]=99;const section=failure(data,'UnsupportedSectionVersion');assert.equal(section.detail.section,1);assert.equal(section.detail.version,99);
 data=bytes.slice();data[19]=255;assert.equal(failure(data,'InvalidArchive').message,'archive 2 block codec mismatch');
 data=bytes.slice();data[18]=0;assert.equal(failure(data,'InvalidArchive').message,'archive 2 required block flag absent');
 data=bytes.slice();data[14]=99;assert.equal(failure(data,'UnknownRequiredSection').detail.section,99);
 data=bytes.slice(0,14);data[12]=0;data[13]=0;assert.equal(failure(data,'MissingRequiredSection').detail.section,1);
 const view=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength),second=24+view.getUint32(20,true);data=bytes.slice();data[second]=1;data[second+1]=0;assert.equal(failure(data,'DuplicateSection').detail.section,1);
 assert.deepEqual(m.toBinary(),bytes);
});
test('Binary input requires Uint8Array and respects subarray boundaries',()=>{
 for(const input of [undefined,null,[],[1,2],new ArrayBuffer(4),new Uint16Array(4),new Int8Array(4),new DataView(new ArrayBuffer(4)),{},'bytes'])assert.throws(()=>b.Molecule.fromBinary(input),TypeError);
 const source=b.Molecule.fromSmiles('CO'),bytes=source.toBinary(),backing=new Uint8Array(bytes.length+8);backing.set(bytes,3);const sub=backing.subarray(3,3+bytes.length),m=b.Molecule.fromBinary(sub);assert.equal(m.toSmiles(),source.toSmiles());assert.deepEqual(m.toBinary(),bytes);assert.equal(typeof Object.getOwnPropertyDescriptor(b.PickleError.prototype,'kind').get,'function');
});
