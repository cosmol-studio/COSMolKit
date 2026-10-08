import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('Four molecule kekulization methods retain all policies, value semantics and reusable parameters',()=>{
 const defaults=new b.KekulizeParams();
 assert.deepEqual([defaults.markAtomsBonds,defaults.canonical,defaults.maxBacktracks],[true,true,100]);
 for(const text of ['c1ccccc1','c1ccncc1','c1ccc2ccccc2c1','CCO']){
  const original=b.Molecule.fromSmiles(text),before=original.toSmiles();
  const value=original.withKekulizedBonds(),explicit=original.withKekulizedBondsWithParams(defaults);
  assert.equal(value.toSmiles(),explicit.toSmiles());
  const mutable=b.Molecule.fromSmiles(text);assert.equal(mutable.kekulizeBonds(),undefined);assert.equal(mutable.toSmiles(),value.toSmiles());
  for(const mark of [false,true])for(const canonical of [false,true])for(const max of [0,1,100,4294967295]){
   const p=new b.KekulizeParams(mark,canonical,max),out=original.withKekulizedBondsWithParams(p),m=b.Molecule.fromSmiles(text);
   assert.equal(m.kekulizeBondsWithParams(p),undefined);assert.equal(m.toSmiles(),out.toSmiles());
   assert.equal(original.withKekulizedBondsWithParams(p).toSmiles(),out.toSmiles());
   assert.deepEqual([p.markAtomsBonds,p.canonical,p.maxBacktracks],[mark,canonical,max]);
   assert.equal(original.toSmiles(),before);
  }
  const detached=original.withKekulizedBonds();original.free();assert.equal(detached.toSmiles(),value.toSmiles());
 }
 const benzene=b.Molecule.fromSmiles('c1ccccc1');assert.ok(!benzene.withKekulizedBonds().toSmiles().includes('c'));assert.equal(benzene.toSmiles(),'c1ccccc1');
});
test('Source kekulization failures retain typed payloads and leave input unchanged for all four calls',()=>{
 for(const [text,kind,atoms] of [['c','AromaticAtomOutsideRing',null],['c1cccc1','NotKekulizable',[0,1,2,3,4]]]){
  const m=b.Molecule.fromSmilesWithParams(text,new b.SmilesParseParams(false)),before=m.toSmiles(),params=new b.KekulizeParams();
  for(const call of [()=>m.withKekulizedBonds(),()=>m.withKekulizedBondsWithParams(params),()=>m.kekulizeBonds(),()=>m.kekulizeBondsWithParams(params)]){
   assert.throws(call,e=>{
    assert.equal(e.name,'OperationError');assert.equal(e.kind,'Kekulize');
    const cause=e.cause;assert.equal(cause.name,'KekulizeError');assert.equal(cause.kind,kind);assert.equal(cause.domain,'kekulize');assert.ok(cause.detail instanceof b.KekulizeError);
    const detail=cause.detail;assert.equal(detail.kind,kind);assert.equal(detail.domain,'kekulize');assert.equal(typeof detail.message,'string');assert.equal(detail.cause,null);
    assert.equal(detail.atom,atoms===null?0:null);assert.deepEqual(detail.problemAtoms,atoms);
    for(const key of ['expected','actual','atomCount','field','begin','end','questions','bitWidth','before','after','bond','reason'])assert.equal(detail[key],null);
    if(atoms!==null){const copy=detail.problemAtoms;copy[0]=99;assert.deepEqual(detail.problemAtoms,atoms);}
    return true;
   });assert.equal(m.toSmiles(),before);
  }
 }
});
test('Kekulization rejects invalid host parameter classes, boolean values and integer widths',()=>{
 const m=b.Molecule.fromSmiles('c1ccccc1'),before=m.toSmiles();
 for(const p of [{},null,undefined,new b.RemoveHsParams()])for(const name of ['withKekulizedBondsWithParams','kekulizeBondsWithParams'])assert.throws(()=>m[name](p),TypeError);
 for(const n of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>new b.KekulizeParams(true,true,n),RangeError);
 assert.throws(()=>new b.KekulizeParams('true'),TypeError);assert.throws(()=>new b.KekulizeParams(true,1),TypeError);assert.equal(m.toSmiles(),before);
});
