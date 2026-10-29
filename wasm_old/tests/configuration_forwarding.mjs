import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Short chemistry and notation calls forward both configuration representations',()=>{
 const mol=ck.Molecule.fromSmiles('CCO'),before=mol.toSmiles();
 const write=new ck.SmilesWriteParams();write.allBondsExplicit=true;
 assert.equal(mol.toSmiles(write),mol.toSmilesWithParams(write));
 assert.equal(mol.toSmiles({allBondsExplicit:true}),mol.toSmilesWithParams(write));
 assert.notEqual(mol.toSmiles(write),before);
 const distance=new ck.DistanceMatrixParams();distance.useBondOrder=true;
 assert.deepEqual(mol.distanceMatrix(distance).values(),mol.distanceMatrixWithParams(distance).values());
 assert.deepEqual(mol.distanceMatrix({useBondOrder:true}).values(),mol.distanceMatrixWithParams(distance).values());
 const parse=new ck.SmartsParseParams();parse.mergeHs=true;
 assert.equal(ck.parseSmarts('[C][H]',parse).numAtoms(),ck.parseSmartsWithParams('[C][H]',parse).numAtoms());
 assert.equal(ck.QueryGraph.fromSmarts('[C][H]',{mergeHs:true}).numAtoms(),1);
 const params=new ck.MolBlockWriteParams();params.format=ck.SdfFormat.V3000;
 for(const name of ['toMol','toSdf','toSdf2d']){
  const expected=mol[name+'WithParams'](params);
  assert.equal(mol[name](params),expected);assert.equal(mol[name]({format:ck.SdfFormat.V3000}),expected);
  assert.ok(expected.includes('V3000'));
 }
 for(const configuration of [params,{format:ck.SdfFormat.V3000}]){
  assert.throws(()=>mol.toSdf3d(configuration),error=>error.name==='MolecularIoError');
 }
 assert.throws(()=>mol.toSdf3dWithParams(params),error=>error.name==='MolecularIoError');
 const geometry=mol.with3dConformer({randomSeed:42});
 assert.equal(geometry.toSdf3d(params),geometry.toSdf3dWithParams(params));
 assert.equal(geometry.toSdf3d({format:ck.SdfFormat.V3000}),geometry.toSdf3dWithParams(params));
 assert.equal(mol.toSmiles(),before);
});

test('Short fingerprint calls use supplied parameters instead of silently using defaults',()=>{
 const mol=ck.Molecule.fromSmiles('CCOC(=O)C');
 for(const [prefix,Class,Generator] of [['fingerprintMorgan',ck.MorganFingerprintParams,ck.MorganGeneratorParams],['fingerprintAtomPair',ck.AtomPairFingerprintParams,ck.AtomPairGeneratorParams],['fingerprintTopologicalTorsion',ck.TopologicalTorsionFingerprintParams,ck.TopologicalTorsionGeneratorParams]]){
  const params=new Class(),generator=new Generator();generator.fpSize=128;params.generator=generator;
  assert.equal(mol[prefix](params).nBits(),128);
  const expected=mol[prefix+'WithParams'](params,null).onBits();
  assert.deepEqual(mol[prefix](params).onBits(),expected);
  assert.deepEqual(mol[prefix]({generator}).onBits(),expected);
  assert.throws(()=>mol[prefix]({unknown:true}),TypeError);
  assert.throws(()=>mol[prefix]({generator:'128'}),TypeError);
 }
});

test('Batch short calls preserve configured parsing, query callbacks and error positions',()=>{
 const parse=new ck.SmilesParseParams(),batchParams=new ck.BatchParams(ck.BatchErrorMode.KeepErrors);
 for(const value of [ck.MoleculeBatch.fromSmilesList(['CCO','[','O'],parse,batchParams),ck.MoleculeBatch.fromSmilesList(['CCO','[','O'],{sanitize:true},{errors:ck.BatchErrorMode.KeepErrors})]){
  assert.deepEqual(value.validMask(),[true,false,true]);
  assert.deepEqual(value.errors().map(e=>e.index()),[1]);
  const write=new ck.SmilesWriteParams();write.allBondsExplicit=true;
  let ticks=0;const callback=()=>{ticks++;};
  assert.deepEqual(value.toSmiles(write,{progressCallback:callback}),['C-C-O',null,'O']);assert.equal(ticks,3);
  assert.deepEqual(value.toSmiles({allBondsExplicit:true},new ck.BatchQueryParams()),value.toSmilesWithParams(write,new ck.BatchQueryParams()));
  assert.deepEqual(value.withHydrogens({},{}).validMask(),[true,false,true]);
  const expected=value.fingerprintMorganWithParams(new ck.MorganFingerprintParams(),new ck.BatchQueryParams()).map(fp=>fp?.onBits()??null);
  assert.deepEqual(value.fingerprintMorgan({},{}).map(fp=>fp?.onBits()??null),expected);
  const thrown=new Error('progress callback');assert.throws(()=>value.toSmiles({}, {progressCallback:()=>{throw thrown;}}),e=>e===thrown);
 }
});

test('Force-field and reaction short calls honor nondefault options and retain typed failures',()=>{
 const molecule=ck.Molecule.fromSmiles('CCO').withHydrogens().with3dConformer({randomSeed:42});
 for(const prefix of ['uff','mmff']){
  const Upper=prefix==='uff'?'Uff':'Mmff';
  const evalParams=new ck[Upper+'EvaluationParams']();
  assert.equal(molecule[prefix+'EnergyGradient'](evalParams).energy(),molecule[prefix+'EnergyGradientWithParams'](evalParams).energy());
  assert.deepEqual(molecule[prefix+'EnergyGradient']({}).gradient(),molecule[prefix+'EnergyGradientWithParams'](evalParams).gradient());
  const params=new ck[Upper+'OptimizationParams']();params.maxIterations=2;
  const actual=molecule['with'+Upper+'Optimized']({maxIterations:2}),expected=molecule['with'+Upper+'OptimizedWithParams'](params);
  assert.equal(actual.statusCode(),expected.statusCode());
  assert.deepEqual(actual.molecule().coordinates3d(0),expected.molecule().coordinates3d(0));
 }
 const reaction=ck.Reaction.fromSmirks('[C:1]>>[C:1]',{});
 const options=new ck.ReactionWriteParams();
 assert.equal(reaction.toSmirks(options),reaction.toSmirksWithParams(options));
 assert.equal(reaction.toSmirks({}),reaction.toSmirksWithParams(options));
 assert.throws(()=>ck.parseSmirks('[C:1]>>[C:1]',{unknown:true}),TypeError);
});
