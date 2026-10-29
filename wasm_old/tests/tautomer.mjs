import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Cross-molecule selection and owned callback snapshots preserve canonical scoring',()=>{
 const enol=ck.Molecule.fromSmiles('C=C(C)O'),keto=ck.Molecule.fromSmiles('CC(C)=O');
 const before=[enol.toSmiles(),keto.toSmiles()],params=new ck.TautomerParams();
 assert.equal(ck.canonicalTautomerFromMolecules([enol,keto,enol]).toSmiles(),'CC(C)=O');
 assert.equal(ck.canonicalTautomerFromMoleculesWithParams([enol,keto],params).toSmiles(),'CC(C)=O');
 assert.equal(ck.canonicalTautomerFromMolecules([enol,keto],{scorer:v=>v.toSmiles()==='C=C(C)O'?100:-100}).toSmiles(),'C=C(C)O');
 assert.deepEqual([enol.toSmiles(),keto.toSmiles()],before);
 assert.throws(()=>ck.canonicalTautomerFromMolecules([]),e=>e.name==='OperationError'&&e.cause.name==='TautomerRunError');
 assert.throws(()=>ck.canonicalTautomerFromMolecules([{}]),TypeError);
 let source,progress;
 keto.enumerateTautomers({callback:(s,p)=>{source=s.toOwned();progress=p.toOwned();s.free();p.free();return false;}});
 enol.free();keto.free();params.free();
 assert.equal(source.numAtoms(),4);assert.ok(progress.len()>0);assert.equal(progress.entries()[0][1].numAtoms(),4);
});

test('Tautomer enumeration retains source keys, COW molecules and bounded status',()=>{
 const mol=ck.Molecule.fromSmiles('CC(=O)C'),before=mol.toSmiles();
 const params=new ck.TautomerParams();params.maxTautomers=100;params.setMaxTransforms(200);
 const enumeration=mol.enumerateTautomers(params);
 assert.equal(enumeration.status,ck.TautomerEnumerationStatus.Completed);
 assert.equal(enumeration.len(),2);assert.equal(enumeration.isEmpty(),false);
 assert.deepEqual(enumeration.canonicalSmiles(),['C=C(C)O','CC(C)=O']);
 assert.deepEqual(enumeration.entries().map(([key,value])=>[key,value.toSmiles()]),[['C=C(C)O','C=C(C)O'],['CC(C)=O','CC(C)=O']]);
 assert.equal(enumeration.canonicalTautomer().toSmiles(),'CC(C)=O');
 assert.equal(mol.canonicalTautomerWithParams(params).toSmiles(),'CC(C)=O');
 assert.equal(mol.enumerateTautomers({maxTautomers:1}).status,ck.TautomerEnumerationStatus.MaxTautomersReached);
 assert.equal(enumeration.get(99),null);assert.equal(enumeration.iter().length,2);
 const child=enumeration.get(0);enumeration.free();assert.equal(child.toSmiles(),'C=C(C)O');
 assert.equal(mol.toSmiles(),before);
 const legacy=ck.TautomerParams.v1(),count=legacy.transformCount;
 legacy.maxTautomers=3;legacy.removeSp3Stereo=false;legacy.scoreParams={terms:[]};
 assert.equal(legacy.transformCount,count);assert.equal(legacy.withMaxTautomers(2).transformCount,count);
 assert.equal(legacy.maxTautomers,3);assert.equal(legacy.withMaxTautomers(2).maxTautomers,2);
 assert.throws(()=>{legacy.maxTautomers=-1;},RangeError);assert.equal(legacy.maxTautomers,3);
});

test('Tautomer scores use typed terms and nested configuration without changing source',()=>{
 const mol=ck.Molecule.fromSmiles('CC(=O)C'),before=mol.toSmiles();
 const term=new ck.TautomerScoreTerm('carbonyl','[C]=[O]',7);
 assert.equal(term.name,'carbonyl');assert.equal(term.score,7);
 const params=new ck.TautomerScoreParams([term]);
 assert.equal(mol.tautomerScoreWithParams(params).substructure,7);
 const score=mol.tautomerScore({terms:[term]});
 assert.equal(score.total,score.ring+score.substructure+score.heteroHydrogen);
 params.terms=[];assert.equal(mol.tautomerScore(params).substructure,0);
 assert.ok(ck.defaultTautomerScoreTerms().length>0);assert.equal(mol.toSmiles(),before);
});

test('Tautomer callbacks and scorers retain restricted snapshots and original exceptions',()=>{
 const mol=ck.Molecule.fromSmiles('CC(=O)C');let calls=0,source,progress,entry;
 const callback=(s,p)=>{calls++;source=s;progress=p;entry=p.entries()[0][1];s.tautomerScore();entry.tautomerScore();return false;};
 const params=new ck.TautomerParams();params.callback=callback;
 assert.equal(params.callback,callback);
 const result=mol.enumerateTautomers(params);
 assert.equal(result.status,ck.TautomerEnumerationStatus.Canceled);assert.equal(calls,1);
 params.free();result.free();assert.equal(source.numAtoms(),4);assert.ok(progress.len()>=1);
 assert.equal(entry.atoms().length,4);assert.equal(entry.atom(99),null);
 assert.equal(entry.atomDegree(1),3);assert.equal(entry.atomMetadata().length,4);
 const scorer=v=>{v.tautomerScore();return v.toSmiles()==='C=C(C)O'?100:-100;};
 assert.equal(mol.canonicalTautomer({scorer}).toSmiles(),'C=C(C)O');
 const error=new Error('original callback failure');
 assert.throws(()=>mol.enumerateTautomers({callback:()=>{throw error;}}),e=>e===error);
 assert.throws(()=>mol.canonicalTautomer({scorer:()=>0.5}),RangeError);
 assert.throws(()=>mol.enumerateTautomers({callback:()=>1}),TypeError);
 assert.equal(mol.toSmiles(),'CC(C)=O');
});

test('Tautomer catalogs retain structured native errors and validate host rows',()=>{
 const empty=ck.TautomerParams.fromTransformData([]);assert.equal(empty.transformCount,0);
 assert.throws(()=>ck.TautomerParams.fromTransformData([['bad']]),TypeError);
 assert.throws(()=>ck.TautomerParams.fromTransformData([['bad','[','', '']]),e=>e.name==='TautomerCatalogError'&&e.kind==='Transform'&&e.detail instanceof ck.TautomerCatalogError&&e.cause instanceof Error);
 assert.throws(()=>ck.TautomerParams.fromTransformFile('/nonexistent/cosmolkit-tautomer.txt'),e=>e.name==='TautomerCatalogError'&&e.kind==='BadInputFile'&&e.path.endsWith('cosmolkit-tautomer.txt'));
});
