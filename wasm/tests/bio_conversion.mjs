import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const pdb="ATOM      1  CA AALA A   1       1.000   2.000   3.000  0.60 20.00           C  \nATOM      2  CA BALA A   1       4.000   5.000   6.000  0.40 30.00           C  \nATOM      3  N   ALA A   1       7.000   8.000   9.000  1.00 10.00           N  \nEND\n",hydrogen="ATOM      1  CA AALA A   1       1.000   2.000   3.000  0.60 20.00           C  \nATOM      2  CA BALA A   1       4.000   5.000   6.000  0.40 30.00           C  \nATOM      3  N   ALA A   1       7.000   8.000   9.000  1.00 10.00           N  \nATOM      4  H   ALA A   1       7.800   8.000   9.000  1.00 20.00           H  \nEND\n",invalid="ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nATOM      2  F   ALA A   1       1.300   0.000   0.000  1.00 20.00           F  \nATOM      3  F   ALA A   1      -1.300   0.000   0.000  1.00 20.00           F  \nATOM      4  F   ALA A   1       0.000   1.300   0.000  1.00 20.00           F  \nATOM      5  F   ALA A   1       0.000  -1.300   0.000  1.00 20.00           F  \nATOM      6  F   ALA A   1       0.000   0.000   1.300  1.00 20.00           F  \nCONECT    1    2    3    4    5    6\nEND\n";
test('BIO full conversion parameters retain defaults, strict widths and reusable class values',()=>{
 const p=new b.BioMoleculeParams();assert.deepEqual([p.sanitize,p.removeHs,p.flavor,p.proximityBonding],[true,true,0,true]);
 const q=new b.BioMoleculeParams(false,false,4294967295,false);assert.deepEqual([q.sanitize,q.removeHs,q.flavor,q.proximityBonding],[false,false,4294967295,false]);
 for(const value of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>new b.BioMoleculeParams(true,true,value),RangeError);
 assert.throws(()=>new b.BioMoleculeParams(true,true,'1'),TypeError);assert.throws(()=>new b.BioMoleculeParams(1),TypeError);assert.throws(()=>new b.BioMoleculeParams(true,null),TypeError);assert.throws(()=>new b.BioMoleculeParams(true,true,0,1),TypeError);
});
test('BIO and Protein produce detached Molecules with alternate locations, hydrogen policy and exact coordinate order',()=>{
 for(const T of [b.BioStructure,b.Protein]){
  const s=T.fromPdb(pdb),normal=s.toMolecule(),params=new b.BioMoleculeParams(true,true,1,false);
  assert.equal(normal.numAtoms(),2);assert.deepEqual(Array.from(normal.coordinates3d(0)),[1,2,3,7,8,9]);
  for(let i=0;i<2;i++){const all=s.toMoleculeWithParams(params);assert.equal(all.numAtoms(),3);assert.equal(all.numBonds(),0);assert.deepEqual(Array.from(all.coordinates3d(0)),[1,2,3,4,5,6,7,8,9]);assert.equal(params.flavor,1);}
  assert.equal(Array.from(s.atoms()).length,3);for(const value of [{},null,undefined,new b.BioReadParams()])assert.throws(()=>s.toMoleculeWithParams(value),TypeError);
  s.free();assert.equal(normal.numAtoms(),2);
  // Conversion returns a detached molecule, not a prepared runtime cache.
  assert.throws(()=>normal.addHydrogens(),e=>e.name==='OperationError'&&e.kind==='Hydrogen'&&e.cause.name==='HydrogenError'&&e.cause.kind==='Valence');
  assert.equal(normal.numAtoms(),2);assert.deepEqual(Array.from(normal.coordinates3d(0)),[1,2,3,7,8,9]);
  const prepared=normal.withAssignedValence();prepared.addHydrogens();assert.equal(prepared.numAtoms()>2,true);assert.equal(normal.numAtoms(),2);
  const h=T.fromPdb(hydrogen),withHs=h.toMoleculeWithParams(new b.BioMoleculeParams(true,false)),withoutHs=h.toMolecule();assert.equal(withHs.numAtoms(),3);assert.equal(withoutHs.numAtoms(),2);assert.equal(Array.from(h.atoms()).length,4);
  const unsanitized=h.toMoleculeWithParams(new b.BioMoleculeParams(false,true));assert.equal(unsanitized.numAtoms(),3);
  assert.equal(T.fromPdb('END\n').toMolecule().numAtoms(),0);
 }
});
test('BIO failed conversion preserves typed owner causes and permits explicitly unsanitized conversion',()=>{
 for(const T of [b.BioStructure,b.Protein]){
  const s=T.fromPdb(invalid);let captured;
  assert.throws(()=>s.toMolecule(),e=>{captured=e;return e instanceof Error&&e.name==='BioMoleculeError'&&e.domain==='bio'&&e.kind==='Conversion'&&e.detail instanceof b.BioMoleculeError&&e.cause instanceof Error&&e.cause.detail instanceof b.BioMoleculeConversionError;});
  assert.equal(captured.detail.cause,captured.cause);assert.equal(captured.cause.domain,'bio');assert.equal(captured.cause.detail.kind,captured.cause.kind);assert.equal(captured.cause.detail.message,captured.cause.message);assert.equal(captured.cause.kind,'Hydrogens');assert.match(captured.cause.message,/Explicit valence for atom # 0 C, 5/);assert.equal(captured.cause.cause,undefined);assert.equal(captured.cause.detail.cause,null);
  assert.equal(Array.from(s.atoms()).length,6);assert.equal(s.toMoleculeWithParams(new b.BioMoleculeParams(false)).numAtoms(),6);
 }
});
