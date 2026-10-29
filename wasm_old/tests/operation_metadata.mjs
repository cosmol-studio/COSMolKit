import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Concrete reader errors and property text remain typed',()=>{
 assert.throws(()=>ck.Molecule.fromSdf(''),e=>e.cause.name==='SdfReadError'&&e.cause.kind==='Empty'&&e.cause.detail instanceof ck.SdfReadError);
 const mol=ck.Molecule.fromSmiles('C').toBuilder().withProperty('note','retained').build();
 assert.equal(ck.propertyValueToText(mol.property('note')),'retained');assert.throws(()=>ck.propertyValueToText({}),TypeError);
 assert.equal(typeof ck.SmilesStereoError,'function');
});

test('Generated declarations remain read-only and all metadata lookups agree',()=>{
 const iterator=ck.featureSpecs(),features=[];
 for(let row;(row=iterator.next())!==null;)features.push(row);
 assert.equal(iterator.next(),null);assert.ok(features.length>0);
 assert.equal(new Set(features.map(v=>v.name)).size,features.length);
 for(const feature of features){
  assert.ok(feature instanceof ck.FeatureSpec);
  const found=ck.featureSpec(feature.name);
  assert.equal(found.category,feature.category);assert.equal(found.docs,feature.docs);
  assert.throws(()=>{feature.name='changed';},TypeError);
 }
 assert.equal(ck.featureSpec('not-a-feature'),null);
 const operations=ck.operationSpecs();assert.ok(operations.length>0);
 for(const row of operations){
  const found=ck.operationSpec(row.method);
  assert.equal(found.implFn,row.implFn);assert.equal(found.resultType,row.resultType);
  assert.equal(found.output,row.output);assert.equal(found.domain,row.domain);
  assert.ok(row.status instanceof ck.FunctionStatus);
  assert.ok(['Parity','ParityWithDifferences','Native','Experimental'].includes(row.status.kind));
  assert.ok(Object.values(ck.ParityPolicy).includes(row.parity));
  assert.equal(typeof row.ioRoundtrip,'boolean');
  assert.equal(typeof row.semanticPreconditions,'number');
  const access=row.access;assert.equal(typeof access.read,'number');assert.equal(typeof access.write,'number');
  access.write=255;assert.deepEqual(row.access,found.access);
  const effects=row.derivedEffects;assert.equal(typeof effects.operationDefined,'number');
  effects.recompute=65535;assert.deepEqual(row.derivedEffects,found.derivedEffects);
  assert.throws(()=>{row.mayMutate=255;},TypeError);
 }
 for(const row of ck.supportMatrix()){
  assert.ok(ck.featureSpec(row.feature.name) instanceof ck.FeatureSpec);
  if(row.operation!==null)assert.equal(ck.operationSpec(row.operation.method).method,row.operation.method);
 }
 for(const row of ck.operationInvariantMatrix())assert.equal(ck.operationInvariant(row.operation.method).profile,row.profile);
 for(const row of ck.parityMatrix()){
  const found=ck.operationParity(row.operation.method);
  assert.equal(found.profile,row.profile);assert.equal(found.rdkitVersion,row.rdkitVersion);
  assert.equal(found.feature.name,row.feature.name);
 }
 for(const lookup of [ck.operationSpec,ck.operationInvariant,ck.operationParity])assert.equal(lookup('not-an-operation'),null);
 iterator.free();assert.equal(ck.featureSpec(features[0].name).name,features[0].name);
});
