import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('typed atom metadata preserves native values and COW',()=>{
 const source=b.Molecule.fromSmiles('CCO');
 for(const [suffix,value] of [['Int',42],['Uint',4294967295],['Bool',true],['Bool',false],['Double',1.25],['Double',-0],['String','α\u0000β'],['IntVector',[1,-2]],['StringVector',['a','β']],['IntVector',[]],['StringVector',[]]]){
  const tagged=source[`withAtomProperty${suffix}`](0,'tracking_id',value);
  assert.deepEqual(tagged[`atomProperty${suffix}`](0,'tracking_id'),value);
  assert.equal(source[`atomProperty${suffix}`](0,'tracking_id'),null);
  const copied=tagged.withAtomPropertyString(0,'copy_note','copy'); tagged.setAtomPropertyInt_(0,'tracking_id',7);
  assert.deepEqual(copied[`atomProperty${suffix}`](0,'tracking_id'),value);
  assert.equal(tagged.atomPropertyInt(0,'tracking_id'),7);
 }
 assert.throws(()=>source.setAtomPropertyString_(0,'_CIPCode','R'),e=>e.kind==='ReservedAtomPropertyKey');
 assert.throws(()=>source.setAtomPropertyInt_(9,'tracking_id',1),e=>e.kind==='AtomPropertyIndex');
 assert.throws(()=>source.setAtomPropertyInt_(0,'tracking_id',{}),TypeError);
});
const sdf="binding record\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n>  <NOTE>\nfirst\nline\n\n>  <NOTE>\nsecond\n\n>  <atom.iprop.score>\n-2147483648 2147483647\n\n>  <atom.prop.label>\nalpha n/a\n\n>  <atom.dprop.weight>\n0.1 -0\n\n>  <atom.bprop.flag>\n1 0\n\n>  <bond.iprop.edge>\n7\n\n>  <bond.prop.label>\nbeta\n\n>  <bond.dprop.weight>\n0.1\n\n>  <bond.bprop.flag>\n1\n\n$$$$\n";
test('Atom and bond property strings use canonical floating/integer/bool spelling and exact null absence',()=>{
 const m=b.Molecule.fromSdf(sdf),before=Array.from(m.coordinates2d());
 for(let repeat=0;repeat<2;repeat++){
  assert.equal(m.atomPropertyInt(0,'score'),-2147483648);assert.equal(m.atomPropertyInt(1,'score'),2147483647);assert.equal(m.atomPropertyString(0,'label'),'alpha');assert.equal(m.atomPropertyString(1,'label'),null);assert.equal(m.atomPropertyDouble(0,'weight'),0.1);assert.equal(m.atomPropertyDouble(1,'weight'),-0);assert.equal(m.atomPropertyBool(0,'flag'),true);assert.equal(m.atomPropertyBool(1,'flag'),false);
  assert.equal(m.bondPropertyString(0,'edge'),'7');assert.equal(m.bondPropertyString(0,'label'),'beta');assert.equal(m.bondPropertyString(0,'weight'),'0.10000000000000001');assert.equal(m.bondPropertyString(0,'flag'),'1');
  for(const id of [0,1,4294967295]){if(id<2)assert.equal(m.atomPropertyString(id,'missing'),null);else assert.throws(()=>m.atomPropertyString(id,'missing'),e=>e.cause.kind==='AtomPropertyIndex');assert.equal(m.bondPropertyString(id,'missing'),null);}assert.throws(()=>m.atomPropertyString(0,'score'),e=>e.variant==='Type');assert.equal(m.bondPropertyString(4294967295,'edge'),null);
 }
 assert.equal(m.numAtoms(),2);assert.equal(m.numBonds(),1);assert.deepEqual(Array.from(m.coordinates2d()),before);
});
test('Property string indices preserve full wasm32 usize and reject host coercion',()=>{
 const m=b.Molecule.fromSdf(sdf);for(const name of ['atomPropertyString','bondPropertyString']){for(const value of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>m[name](value,'score'),RangeError);for(const value of ['0',0n,null,{}])assert.throws(()=>m[name](value,'score'),TypeError);if(name==='bondPropertyString')assert.equal(m[name](4294967295,'missing'),null);else assert.throws(()=>m[name](4294967295,'missing'));}
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
