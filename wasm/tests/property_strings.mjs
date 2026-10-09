import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('typed atom metadata preserves native values and COW',()=>{
 const source=b.Molecule.fromSmiles('CCO');
 for(const value of [42,4294967295,true,false,1.25,-0,'α\u0000β',[1,-2],['a','β'],[]]){
  const tagged=source.withAtomProperty(0,'tracking_id',value);
  assert.deepEqual(tagged.atomProperty(0,'tracking_id'),value);
  assert.equal(source.atomProperty(0,'tracking_id'),null);
  const copied=tagged.withAtomProperty(0,'copy_note','copy'); tagged.setAtomProperty(0,'tracking_id',7);
  assert.deepEqual(copied.atomProperty(0,'tracking_id'),value);
  assert.equal(tagged.atomProperty(0,'tracking_id'),7);
 }
 assert.throws(()=>source.setAtomProperty(0,'_CIPCode','R'),e=>e.kind==='ReservedAtomPropertyKey');
 assert.throws(()=>source.setAtomProperty(9,'tracking_id',1),e=>e.kind==='AtomPropertyIndex');
 assert.throws(()=>source.setAtomProperty(0,'tracking_id',{}),TypeError);
});
const sdf="binding record\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n>  <NOTE>\nfirst\nline\n\n>  <NOTE>\nsecond\n\n>  <atom.iprop.score>\n-2147483648 2147483647\n\n>  <atom.prop.label>\nalpha n/a\n\n>  <atom.dprop.weight>\n0.1 -0\n\n>  <atom.bprop.flag>\n1 0\n\n>  <bond.iprop.edge>\n7\n\n>  <bond.prop.label>\nbeta\n\n>  <bond.dprop.weight>\n0.1\n\n>  <bond.bprop.flag>\n1\n\n$$$$\n";
test('Atom and bond property strings use canonical floating/integer/bool spelling and exact null absence',()=>{
 const m=b.Molecule.fromSdf(sdf),before=Array.from(m.coordinates2d());
 for(let repeat=0;repeat<2;repeat++){
  assert.equal(m.atomPropertyString(0,'score'),'-2147483648');assert.equal(m.atomPropertyString(1,'score'),'2147483647');assert.equal(m.atomPropertyString(0,'label'),'alpha');assert.equal(m.atomPropertyString(1,'label'),null);assert.equal(m.atomPropertyString(0,'weight'),'0.10000000000000001');assert.equal(m.atomPropertyString(1,'weight'),'-0');assert.equal(m.atomPropertyString(0,'flag'),'1');assert.equal(m.atomPropertyString(1,'flag'),'0');
  assert.equal(m.bondPropertyString(0,'edge'),'7');assert.equal(m.bondPropertyString(0,'label'),'beta');assert.equal(m.bondPropertyString(0,'weight'),'0.10000000000000001');assert.equal(m.bondPropertyString(0,'flag'),'1');
  for(const id of [0,1,4294967295]){assert.equal(m.atomPropertyString(id,'missing'),null);assert.equal(m.bondPropertyString(id,'missing'),null);}assert.equal(m.atomPropertyString(4294967295,'score'),null);assert.equal(m.bondPropertyString(4294967295,'edge'),null);
 }
 assert.equal(m.numAtoms(),2);assert.equal(m.numBonds(),1);assert.deepEqual(Array.from(m.coordinates2d()),before);
});
test('Property string indices preserve full wasm32 usize and reject host coercion',()=>{
 const m=b.Molecule.fromSdf(sdf);for(const name of ['atomPropertyString','bondPropertyString']){for(const value of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>m[name](value,'score'),RangeError);for(const value of ['0',0n,null,{}])assert.throws(()=>m[name](value,'score'),TypeError);assert.equal(m[name](4294967295,'missing'),null);}
});
test('PropertyText projection preserves Unicode and embedded NUL without mutating the receiver',()=>{
 const source=b.Molecule.fromSmiles('CCO');
 const value=source.withName('乙醇\u0000tail').withProperty('label','α\u0000β');
 assert.equal(value.nameOrEmpty(),'乙醇\u0000tail');
 assert.equal(value.propertyOrEmpty('label'),'α\u0000β');
 assert.equal(source.nameOrEmpty(),'');
 assert.equal(source.propertyOrEmpty('label'),'');
 const properties=b.SdfRecord.fromSdf(sdf).properties();
 const names=properties.prop('__computedProps');
 assert.equal(names.kind(),b.PropertyValueKind.StringVector);
 assert.deepEqual(names.asStringVector(),properties.computedPropNames());
 assert.throws(()=>names.asInt(),e=>e.kind==='KindMismatch'&&e.actual===b.PropertyValueKind.StringVector);
});
