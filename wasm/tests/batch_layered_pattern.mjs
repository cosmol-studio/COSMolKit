import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const batch=()=>b.MoleculeBatch.fromSmilesListWithParams(['CCO','[','O'],new b.SmilesParseParams(),new b.BatchParams(b.BatchErrorMode.KeepErrors));
test('Layered and Pattern batch outputs equal scalar outputs', () => {
 for (const text of ['CCO', 'c1ccncc1', 'CC(=O)N']) {
  const molecule = b.Molecule.fromSmiles(text), values = b.MoleculeBatch.fromSmilesList([text]);
  assert.deepEqual(molecule.fingerprintLayered().onBits(), values.fingerprintLayered()[0].onBits());
  assert.deepEqual(molecule.fingerprintPattern().onBits(), values.fingerprintPattern()[0].onBits());
  molecule.free(); values.free();
 }
});
test('Layered/Pattern parameters retain every default, unknown flags, borrowed masks and nullable roots',()=>{
 const p=new b.LayeredFingerprintParams();assert.deepEqual([p.layers,p.minPath,p.maxPath,p.fpSize,p.atomCounts,p.setOnlyBits,p.branchedPaths,p.fromAtoms],[4294967295,1,7,2048,null,null,true,null]);
 assert.equal(b.LayeredFingerprintLayers.fromBitsRetain(4294967295).bits(),4294967295);assert.equal(b.LayeredFingerprintLayers.fromBitsRetain(2147483648).bits(),2147483648);
 const mask=b.Fingerprint.fromOnBits(2048,[674]),seed=[10,20,30],roots=[0],custom=new b.LayeredFingerprintParams(63,2,5,2048,seed,mask,false,roots);seed[0]=99;roots.push(1);
 assert.deepEqual([custom.layers,custom.minPath,custom.maxPath,custom.fpSize,custom.branchedPaths],[63,2,5,2048,false]);assert.deepEqual(custom.atomCounts,[10,20,30]);assert.deepEqual(custom.fromAtoms,[0]);assert.deepEqual(custom.setOnlyBits.onBits(),[674]);assert.deepEqual(mask.onBits(),[674]);mask.free();assert.deepEqual(custom.setOnlyBits.onBits(),[674]);
 assert.deepEqual(new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,[],null,true,[]).fromAtoms,[]);
 const pattern=new b.PatternFingerprintParams();assert.equal(pattern.nBits,2048);assert.equal(pattern.tautomeric,false);assert.deepEqual([new b.PatternFingerprintParams(128,true).nBits,new b.PatternFingerprintParams(128,true).tautomeric],[128,true]);
 for(const n of [-1,4294967296,1.5,NaN]){assert.throws(()=>b.LayeredFingerprintLayers.fromBitsRetain(n),RangeError);assert.throws(()=>new b.PatternFingerprintParams(n),RangeError);}
 assert.throws(()=>new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,null,{}),TypeError);assert.throws(()=>new b.PatternFingerprintParams(undefined,'false'),TypeError);{const previous=p.layers;p.layers=0;assert.deepEqual(p.layers,0);p.layers=previous;}
});
test('all six Layered/Pattern batch methods preserve null rows, detached results and exact source count policy',()=>{
 const value=batch(),query=new b.BatchQueryParams(),options=new b.LayeredFingerprintParams(),pattern=new b.PatternFingerprintParams();
 for(const [name,params,Class,read] of [['fingerprintLayered',options,'Fingerprint',v=>v.onBits()],['fingerprintLayeredWithOutput',options,'LayeredFingerprintResult',v=>v.fingerprint().onBits()],['fingerprintPattern',pattern,'Fingerprint',v=>v.onBits()]]){const implicit=value[name](),explicit=value[name+'WithParams'](params,query);assert.equal(implicit.length,3);assert.equal(implicit[1],null);assert.ok(implicit[0] instanceof b[Class]);assert.deepEqual(implicit.map(v=>v===null?null:read(v)),explicit.map(v=>v===null?null:read(v)));assert.deepEqual(b.MoleculeBatch.fromSmilesList([])[name](),[]);}
 assert.equal(value.fingerprintLayeredWithOutput()[0].atomCounts(),null);
 // Pinned owner layered.rs::original_masks_empty_roots_and_rooted_linear_results_are_preserved.
 const single=b.MoleculeBatch.fromSmilesList(['CCO']),mask=b.Fingerprint.fromOnBits(2048,[674]),seeded=new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,[10,20,30],mask),result=single.fingerprintLayeredWithOutputWithParams(seeded,query)[0];
 assert.deepEqual(result.fingerprint().onBits(),[674]);assert.deepEqual(result.atomCounts(),[11,22,31]);assert.deepEqual(seeded.atomCounts,[10,20,30]);const counts=result.atomCounts();counts[0]=999;assert.deepEqual(result.atomCounts(),[11,22,31]);
 const empty=single.fingerprintLayeredWithOutputWithParams(new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,[0,0,0],null,true,[]),query)[0];assert.deepEqual(empty.fingerprint().onBits(),[]);assert.deepEqual(empty.atomCounts(),[0,0,0]);
 assert.equal(single.fingerprintPatternWithParams(new b.PatternFingerprintParams(128,true),query)[0].nBits(),128);assert.deepEqual(value.validMask(),[true,false,true]);
});
test('Layered/Pattern typed source failures and callback identity remain observable',()=>{
 const value=batch();let ticks=0;const query=new b.BatchQueryParams(null,false,()=>{ticks++;});value.fingerprintLayeredWithParams(new b.LayeredFingerprintParams(),query);assert.equal(ticks,3);ticks=0;value.fingerprintPatternWithParams(new b.PatternFingerprintParams(),query);assert.equal(ticks,3);
 const thrown=Symbol('progress');assert.throws(()=>value.fingerprintLayeredWithParams(new b.LayeredFingerprintParams(),new b.BatchQueryParams(null,false,()=>{throw thrown;})),e=>e===thrown);
 for(const [options,reason] of [[new b.LayeredFingerprintParams(undefined,0),'minPath==0'],[new b.LayeredFingerprintParams(undefined,5,4),'maxPath<minPath'],[new b.LayeredFingerprintParams(undefined,undefined,undefined,0),'fpSize==0']]){ticks=0;assert.throws(()=>value.fingerprintLayeredWithParams(options,query),e=>e.cause.cause.kind==='InvalidArguments'&&e.cause.cause.reason===reason&&e.cause.cause.domain==='Fingerprint'&&e.cause.cause.detail instanceof b.LayeredFingerprintError);assert.equal(ticks,3);}
 assert.throws(()=>value.fingerprintPatternWithParams(new b.PatternFingerprintParams(0),query),e=>e.cause.cause.kind==='EmptyFingerprint'&&e.cause.cause.domain==='Fingerprint'&&e.cause.cause.detail instanceof b.PatternFingerprintError);
 const single=b.MoleculeBatch.fromSmilesList(['CCO']);for(const [options,reason] of [[new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,[0]),'bad atomCounts size'],[new b.LayeredFingerprintParams(undefined,undefined,undefined,undefined,null,b.Fingerprint.fromOnBits(64,[])),'bad setOnlyBits size']])assert.throws(()=>single.fingerprintLayeredWithOutputWithParams(options,new b.BatchQueryParams()),e=>e.recordErrors[0].index()===0&&e.cause.cause.reason===reason);
});
