import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const batch=()=>b.MoleculeBatch.fromSmilesListWithParams(['CCO','[','O'],new b.SmilesParseParams(),new b.BatchParams(b.BatchErrorMode.KeepErrors));
test('AtomPair complete parameter fields preserve defaults, borrowed objects and nullable arrays',()=>{
 const p=new b.AtomPairParams(),fields=['minDistance','maxDistance','includeChirality','use2d','countSimulation','fpSize','bitsPerFeature'];
 assert.deepEqual(fields.map(n=>p[n]),[1,30,false,true,true,2048,1]);assert.deepEqual(p.countBounds,[1,2,4,8]);
 const bounds=[1,3],custom=new b.AtomPairParams(2,20,true,false,false,128,2,bounds);bounds.push(8);assert.deepEqual(custom.countBounds,[1,3]);
 assert.deepEqual(fields.map(n=>custom[n]),[2,20,true,false,false,128,2]);assert.deepEqual(new b.AtomPairParams(undefined,undefined,undefined,undefined,undefined,undefined,undefined,[]).countBounds,[]);
 assert.equal(typeof p.infoString(),'string');assert.equal(p.withJson(p.toJson()).toJson(),p.toJson());assert.throws(()=>p.withJson('{'),e=>e.kind==='Parse'&&e.detail instanceof b.FingerprintJsonError);
 const invariants=new b.AtomPairAtomInvariantsGenerator(true,true);assert.equal(invariants.includeChirality,true);assert.equal(invariants.topologicalTorsionCorrection,true);assert.equal(typeof invariants.infoString(),'string');assert.equal(typeof invariants.toJson(),'string');
 const atoms=[0,2],call=new b.AtomPairFingerprintParams(p,atoms,[],[6,6,8],[],2,invariants,false);atoms.push(1);
 assert.deepEqual(call.fromAtoms,[0,2]);assert.deepEqual(call.ignoreAtoms,[]);assert.deepEqual(call.customAtomInvariants,[6,6,8]);assert.deepEqual(call.customBondInvariants,[]);assert.equal(call.conformerId,2);assert.equal(call.useLegacyStereoPerception,false);assert.equal(call.atomInvariantsGenerator.includeChirality,true);
 assert.equal(p.fpSize,2048);assert.equal(invariants.includeChirality,true);assert.equal(call.generator.fpSize,2048);const defaults=new b.AtomPairFingerprintParams();assert.equal(defaults.fromAtoms,null);assert.equal(defaults.ignoreAtoms,null);assert.equal(defaults.customAtomInvariants,null);assert.equal(defaults.customBondInvariants,null);assert.equal(defaults.atomInvariantsGenerator,null);assert.equal(defaults.conformerId,-1);assert.equal(defaults.useLegacyStereoPerception,true);
 assert.throws(()=>{p.fpSize=9;},TypeError);assert.throws(()=>new b.AtomPairFingerprintParams({}),TypeError);assert.throws(()=>new b.AtomPairParams(-1),RangeError);assert.throws(()=>new b.AtomPairParams(undefined,undefined,'false'),TypeError);
 for(const n of [-2147483649,2147483648,1.5])assert.throws(()=>new b.AtomPairFingerprintParams(null,null,null,null,null,n),RangeError);
});
test('all ten batch AtomPair calls preserve order, real output types and sparse value shapes',()=>{
 const value=batch(),options=new b.AtomPairFingerprintParams(),query=new b.BatchQueryParams();
 const families=[['fingerprintAtomPairList','Fingerprint','onBits'],['fingerprintAtomPairSparseCountList','SparseCountFingerprint','nonzeroElements'],['fingerprintAtomPairCountList','SparseCountFingerprint32','nonzeroElements'],['fingerprintAtomPairSparseBitsList','SparseBitFingerprint','onBits']];
 for(const [name,Class,read] of families){const implicit=value[name](),explicit=value[name+'WithParams'](options,query);assert.equal(implicit.length,3);assert.equal(implicit[1],null);assert.ok(implicit[0] instanceof b[Class]);assert.deepEqual(explicit.map(v=>v?.[read]()??null),implicit.map(v=>v?.[read]()??null));assert.deepEqual(b.MoleculeBatch.fromSmilesList([])[name](),[]);}
 const sparse=value.fingerprintAtomPairSparseCountList()[0];assert.equal(typeof sparse.length(),'bigint');for(const [key,v] of sparse.nonzeroElements()){assert.equal(typeof key,'bigint');assert.equal(typeof v,'number');}
 const signed=value.fingerprintAtomPairSparseBitsList()[0];assert.equal(typeof signed.nBits(),'number');assert.ok(signed.onBits().every(Number.isInteger));
 const outputs=value.fingerprintAtomPairWithOutputList(),explicit=value.fingerprintAtomPairWithOutputListWithParams(options,true,query);assert.equal(outputs[1],null);assert.ok(outputs[0] instanceof b.BatchFingerprintOutput);assert.deepEqual(outputs[0].fingerprint().onBits(),explicit[0].fingerprint().onBits());assert.deepEqual(outputs[0].fingerprint().onBits(),value.fingerprintAtomPairList()[0].onBits());
 assert.deepEqual(value.validMask(),[true,false,true]);
});
test('batch fingerprint additional output is detached and missing collection retains a typed failure',()=>{
 const value=batch(),options=new b.AtomPairFingerprintParams(),query=new b.BatchQueryParams(),result=value.fingerprintAtomPairWithOutputList()[0],out=result.additionalOutput();
 assert.ok(out instanceof b.BatchFingerprintAdditionalOutput);const counts=out.atomCounts();assert.equal(counts.length,3);counts[0]=999;assert.notEqual(out.atomCounts()[0],999);
 const atomBits=out.atomToBits();assert.equal(atomBits.length,3);assert.ok(atomBits.flat().every(v=>typeof v==='bigint'));atomBits[0].push(999n);assert.ok(!out.atomToBits()[0].includes(999n));
 const info=out.bitInfoMap();assert.ok(info instanceof Map);for(const [key,pairs] of info){assert.equal(typeof key,'bigint');assert.ok(pairs.every(v=>v.length===2&&v.every(Number.isInteger)));}info.set(999n,[[9,9]]);assert.ok(!out.bitInfoMap().has(999n));
 assert.equal(out.bitPaths(),null);assert.ok(out.atomsPerBit() instanceof Map);
 const missing=value.fingerprintAtomPairWithOutputListWithParams(options,false,query)[0];assert.throws(()=>missing.additionalOutput(),e=>e.name==='BatchFingerprintOutputError'&&e.kind==='MissingAdditionalOutput'&&e.fingerprintKind==='AtomPair'&&e.detail instanceof b.BatchFingerprintOutputError);
 assert.throws(()=>value.fingerprintAtomPairWithOutputListWithParams(options,'false',query),TypeError);
});
test('AtomPair batch callbacks and source failures preserve row indices and original exception identity',()=>{
 const value=batch(),options=new b.AtomPairFingerprintParams();let ticks=0;const query=new b.BatchQueryParams(null,false,()=>{ticks++;});value.fingerprintAtomPairListWithParams(options,query);assert.equal(ticks,3);
 const thrown=Symbol('callback');ticks=0;assert.throws(()=>value.fingerprintAtomPairListWithParams(options,new b.BatchQueryParams(null,false,()=>{ticks++;throw thrown;})),e=>e===thrown);assert.equal(ticks,3);
 const zero=new b.AtomPairFingerprintParams(new b.AtomPairParams(undefined,undefined,undefined,undefined,undefined,0));ticks=0;assert.throws(()=>value.fingerprintAtomPairListWithParams(zero,query),e=>e.domain==='batch'&&e.cause.cause.kind==='EmptyFingerprint');assert.equal(ticks,0);
 // Both switches must be off: sanitize=false alone retains the default
 // remove-hydrogens path and is not a promise of an uninitialized cache.
 const unprepared=b.Molecule.fromSmilesWithParams('CC',new b.SmilesParseParams(false,undefined,undefined,undefined,false)),raw=b.MoleculeBatch.fromRecords([b.BatchRecord.molecule(unprepared)],b.BatchErrorMode.Strict);
 assert.throws(()=>raw.fingerprintAtomPairList(),e=>{assert.deepEqual(e.recordErrors.map(v=>v.index()),[0]);assert.equal(e.cause.cause.kind,'Preparation');assert.ok(e.cause.cause.detail instanceof b.AtomPairReadError);assert.equal(e.cause.cause.cause.kind,'MissingPreparedValence');assert.ok(e.cause.cause.cause.detail instanceof b.FingerprintPreparationError);return true;});
});
