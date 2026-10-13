import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Potential stereo preserves typed results, clean snapshots and original molecule',()=>{
 const mol=ck.Molecule.fromSmiles('CC(O)Cl'),before=mol.toSmiles();
 const params=new ck.PotentialStereoParams();params.clean=true;params.flagPossible=true;
 const result=mol.potentialStereo(params),explicit=mol.potentialStereoWithParams(params);
 assert.equal(result.stereo.length,1);assert.equal(explicit.stereo.length,1);
 const row=result.stereo[0];
 assert.equal(row.stereoType,ck.PotentialStereoType.AtomTetrahedral);
 assert.equal(row.specified,ck.PotentialStereoSpecified.Unspecified);
 assert.equal(row.centeredOn.kind,'atom');assert.equal(row.centeredOn.index,1);
 assert.equal(row.descriptor,ck.PotentialStereoDescriptor.None);
 // PotentialStereoInfo contains explicit graph neighbors, not the implicit H.
 assert.deepEqual(row.controllingAtoms,[0,2,3]);
 assert.ok(result.cleanedMolecule instanceof ck.Molecule);
 assert.equal(mol.potentialStereo().cleanedMolecule,null);
 assert.equal(mol.potentialStereo({flagPossible:false}).stereo.length,0);
 assert.equal(mol.toSmiles(),before);
 mol.free();result.free();assert.equal(row.centeredOn.index,1);
});

test('CIP configuration forwards selected atoms and propagates concrete errors atomically',()=>{
 const mol=ck.Molecule.fromSmiles('C[C@H](O)Cl'),before=mol.toSmiles();
 const options=new ck.CipLabelOptions([1],null,1000);
 options.maxRecursiveIterations=2000;options.atoms=[1];
 const assigned=mol.withCipLabels(options);
 // Candidate CIPLabeler.cpp sets the molecule marker only for a full selection.
 assert.equal(assigned.cipComputed(),false);
 assert.equal(assigned.atomProperty(1,'_CIPCode'),'R');
 assert.deepEqual(assigned.findChiralCenters(),[[1,'R']]);
 assert.equal(mol.withCipLabels({atoms:[1],maxRecursiveIterations:2000}).cipComputed(),false);
 mol.assignCipLabelsWithOptions_(options);assert.equal(mol.cipComputed(),false);
 assert.equal(mol.atomProperty(1,'_CIPCode'),'R');
 assert.throws(()=>mol.assignCipLabels_({atoms:[99]}),e=>e.name==='OperationError'&&e.cause.name==='CipLabelerError'&&e.cause.detail instanceof ck.CipLabelerError);
 assert.equal(mol.toSmiles(),before);assert.equal(mol.cipComputed(),false);
 assert.equal(mol.atomProperty(1,'_CIPCode'),'R');
 const complete=mol.withCipLabels();assert.equal(complete.cipComputed(),true);complete.free();
 const frames=assigned.tetrahedralStereo();assert.equal(frames.length,1);
 assert.equal(frames[0].center,1);assert.equal(frames[0].ligands.length,4);
 assert.equal(frames[0].ligands.filter(v=>v.isImplicitHydrogen).length,1);
 const ligand=new ck.LigandRef();assert.equal(ligand.atom,null);assert.equal(ligand.isImplicitHydrogen,true);
 assert.equal(new ck.LigandRef(2).atom,2);
 assigned.perceiveStereochemistry();assigned.free();assert.equal(frames[0].center,1);
});

test('Structure tag parameters retain defaults and actual conformer errors',()=>{
 const source=ck.Molecule.fromSmiles('FC(Cl)(Br)I');
 // Canonical structure tagging is a no-op without any conformer.
 assert.equal(source.withChiralTagsFromStructure({conformerId:99}).toSmiles(),source.toSmiles());
 const mol=source.withAdded3dConformer([[1,0,0],[0,0,0],[0,1,0],[0,0,1],[-1,-1,-1]]),before=mol.toSmiles();
 const params=new ck.StructureTagParams();assert.equal(params.conformerId,-1);
 params.replaceExistingTags=false;params.conformerId=99;
 assert.throws(()=>mol.withChiralTagsFromStructure(params),e=>e.name==='OperationError'&&e.cause.name==='StereoError'&&e.cause.detail instanceof ck.StereoError);
 assert.throws(()=>mol.assignChiralTagsFromStructure_({conformerId:99}),e=>e.name==='OperationError');
 assert.equal(mol.toSmiles(),before);
 assert.throws(()=>{params.conformerId=0.5;},RangeError);
 assert.equal(params.conformerId,99);
});
