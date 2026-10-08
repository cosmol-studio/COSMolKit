import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('Six ring-cache and family operations retain all policies, reusable parameters and value semantics',()=>{
 const defaults=new b.RingSearchParams();assert.deepEqual([defaults.includeDativeBonds,defaults.includeHydrogenBonds],[false,false]);
 for(const [text,rings] of [['',0],['CCO',0],['C1CCCCC1',1],['c1ccccc1',1],['c1ccc2ccccc2c1',2],['C1CC2CCC1C2',2],['C1CC1.C1CCC1',2]]){
  const parse=new b.SmilesParseParams(false,undefined,undefined,undefined,false),m=b.Molecule.fromSmilesWithParams(text,parse),before=m.toSmiles();
  assert.throws(()=>m.numRings(),e=>e.kind==='MissingInitializedRings');
  const assigned=m.withAssignedRings();assert.equal(assigned.numRings(),rings);assert.throws(()=>m.numRings(),e=>e.kind==='MissingInitializedRings');
  assert.equal(m.assignRings(),undefined);assert.equal(m.numRings(),rings);assert.equal(m.toSmiles(),before);
  if(text===''){
   for(const fn of [()=>m.withAssignedRingFamilies(),()=>m.assignRingFamilies(),()=>m.withAssignedRingFamiliesWithParams(defaults),()=>m.assignRingFamiliesWithParams(defaults)])assert.throws(fn,e=>e.name==='OperationError'&&e.kind==='Rings'&&e.cause.message==='graph has no nodes');
   for(const d of [false,true])for(const h of [false,true]){const p=new b.RingSearchParams(d,h);assert.throws(()=>m.withAssignedRingFamiliesWithParams(p),e=>e.kind==='Rings');assert.throws(()=>m.assignRingFamiliesWithParams(p),e=>e.kind==='Rings');}
   assert.equal(m.numRings(),0);const out=m.withAssignedRings();m.free();assert.equal(out.numRings(),0);continue;
  }
  const family=m.withAssignedRingFamilies(),explicit=m.withAssignedRingFamiliesWithParams(defaults);assert.equal(family.toSmiles(),explicit.toSmiles());assert.equal(family.numRings(),rings);assert.equal(m.assignRingFamilies(),undefined);assert.equal(m.numRings(),rings);
  for(const dative of [false,true])for(const hydrogen of [false,true]){
   const params=new b.RingSearchParams(dative,hydrogen),out=m.withAssignedRingFamiliesWithParams(params);
   assert.equal(out.toSmiles(),before);assert.equal(m.assignRingFamiliesWithParams(params),undefined);assert.equal(m.numRings(),rings);assert.equal(m.withAssignedRingFamiliesWithParams(params).toSmiles(),before);assert.deepEqual([params.includeDativeBonds,params.includeHydrogenBonds],[dative,hydrogen]);
  }
  const detached=m.withAssignedRings();m.free();assert.equal(detached.numRings(),rings);
 }
});
test('Ring policy host boundaries reject coercion and incorrect parameter classes without changing cache',()=>{
 const m=b.Molecule.fromSmilesWithParams('C1CCCCC1',new b.SmilesParseParams(false,undefined,undefined,undefined,false)),before=m.toSmiles();
 for(const value of [{},undefined,null,new b.KekulizeParams()])for(const name of ['withAssignedRingFamiliesWithParams','assignRingFamiliesWithParams'])assert.throws(()=>m[name](value),TypeError);
 assert.throws(()=>new b.RingSearchParams('true'),TypeError);assert.throws(()=>new b.RingSearchParams(true,1),TypeError);assert.throws(()=>{new b.RingSearchParams().includeHydrogenBonds=true;},TypeError);
 assert.equal(m.toSmiles(),before);assert.throws(()=>m.numRings(),e=>e.kind==='MissingInitializedRings');
});
