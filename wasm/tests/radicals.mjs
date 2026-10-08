import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('Both radical methods preserve molecule value semantics and canonical unsanitized assignment',()=>{
 const p=new b.SmilesParseParams(false);
 for(const text of ['', 'C','[C]','[CH]','[CH2]','[CH3]','[O]','[OH]','[H]','[Fe]','[Fe+3]','[Na+]','[Cl-]','*','C[CH2]','[O][O]']){
  const m=b.Molecule.fromSmilesWithParams(text,p),before=m.toSmiles(),out=m.withAssignedRadicals();
  assert.equal(m.toSmiles(),before);assert.equal(out.numAtoms(),m.numAtoms());assert.equal(m.assignRadicals(),undefined);assert.equal(m.toSmiles(),out.toSmiles());
  const second=out.withAssignedRadicals();assert.equal(second.toSmiles(),out.toSmiles());
  const detached=m.withAssignedRadicals();m.free();assert.equal(detached.toSmiles(),out.toSmiles());
 }
 const carbon=b.Molecule.fromSmilesWithParams('[CH3]',p);assert.match(carbon.withAssignedRadicals().toMol(),/M  RAD/);assert.doesNotMatch(carbon.toMolWithParams(new b.MolBlockWriteParams(undefined,undefined,undefined,false)),/M  RAD/);
});
