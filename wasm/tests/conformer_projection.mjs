import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const calls=[{"name": "with3dConformer", "params": false, "multiple": false, "report": false, "inplace": false}, {"name": "with3dConformerWithParams", "params": true, "multiple": false, "report": false, "inplace": false}, {"name": "embed3dConformer", "params": false, "multiple": false, "report": false, "inplace": true}, {"name": "embed3dConformerWithParams", "params": true, "multiple": false, "report": false, "inplace": true}, {"name": "with3dConformerResult", "params": false, "multiple": false, "report": true, "inplace": false}, {"name": "with3dConformerResultWithParams", "params": true, "multiple": false, "report": true, "inplace": false}, {"name": "embed3dConformerResult", "params": false, "multiple": false, "report": true, "inplace": true}, {"name": "embed3dConformerResultWithParams", "params": true, "multiple": false, "report": true, "inplace": true}, {"name": "with3dConformers", "params": false, "multiple": true, "report": false, "inplace": false}, {"name": "with3dConformersWithParams", "params": true, "multiple": true, "report": false, "inplace": false}, {"name": "embed3dConformers", "params": false, "multiple": true, "report": false, "inplace": true}, {"name": "embed3dConformersWithParams", "params": true, "multiple": true, "report": false, "inplace": true}, {"name": "with3dConformersResult", "params": false, "multiple": true, "report": true, "inplace": false}, {"name": "with3dConformersResultWithParams", "params": true, "multiple": true, "report": true, "inplace": false}, {"name": "embed3dConformersResult", "params": false, "multiple": true, "report": true, "inplace": true}, {"name": "embed3dConformersResultWithParams", "params": true, "multiple": true, "report": true, "inplace": true}];

test('EmbedParams all factories, frozen constructor fields, JSON and checked widths',()=>{
 const p=new b.EmbedParams();assert.equal(p.toJson(),b.EmbedParams.new().toJson());assert.equal(p.toJson(),b.EmbedParams.dg().toJson());assert.equal(p.randomSeed(),-1);assert.equal(p.numThreads(),1);assert.equal(p.maxIterations(),0);assert.equal(p.etVersion(),2);assert.equal(p.clearConfs(),true);assert.equal(p.coordMap(),null);assert.equal(p.cpci(),null);assert.deepEqual([...p.failures()],[]);
 for(const factory of ['new','dg','kdg','etdg','etdgV2','etkdg','etkdgV2','etkdgV3','srEtkdgV3']){const v=b.EmbedParams[factory]();assert.ok(v instanceof b.EmbedParams);assert.equal(v.withJson(v.toJson()).toJson(),v.toJson());}assert.equal(b.EmbedParams.kdg().useBasicKnowledge(),true);assert.equal(b.EmbedParams.etdg().useExpTorsionAnglePrefs(),true);assert.equal(b.EmbedParams.etkdgV3().useMacrocycleTorsions(),true);assert.equal(b.EmbedParams.srEtkdgV3().useSmallRingTorsions(),true);
 const values=[42,-7,-7,!p.clearConfs(),!p.useRandomCoords(),3.25,!p.randNegEig(),42,new Map(),3.25,!p.ignoreSmoothingFailures(),!p.enforceChirality(),!p.useExpTorsionAnglePrefs(),!p.useBasicKnowledge(),!p.verbose(),3.25,3.25,!p.onlyHeavyAtomsForRms(),42,!p.embedFragmentsSeparately(),!p.useSmallRingTorsions(),!p.useMacrocycleTorsions(),!p.useMacrocycle14Config(),42,new Map(),!p.forceTransAmides(),!p.useSymmetryForPruning(),3.25,!p.trackFailures(),!p.enableSequentialRandomSeeds(),!p.symmetrizeConjugatedTerminalGroupsForPruning()];const configured=new b.EmbedParams(...values);
 assert.equal(configured.maxIterations(),values[0]);
 assert.equal(configured.numThreads(),values[1]);
 assert.equal(configured.randomSeed(),values[2]);
 assert.equal(configured.clearConfs(),values[3]);
 assert.equal(configured.useRandomCoords(),values[4]);
 assert.equal(configured.boxSizeMult(),values[5]);
 assert.equal(configured.randNegEig(),values[6]);
 assert.equal(configured.numZeroFail(),values[7]);
 assert.equal(configured.coordMap().size,0);
 assert.equal(configured.optimizerForceTol(),values[9]);
 assert.equal(configured.ignoreSmoothingFailures(),values[10]);
 assert.equal(configured.enforceChirality(),values[11]);
 assert.equal(configured.useExpTorsionAnglePrefs(),values[12]);
 assert.equal(configured.useBasicKnowledge(),values[13]);
 assert.equal(configured.verbose(),values[14]);
 assert.equal(configured.basinThresh(),values[15]);
 assert.equal(configured.pruneRmsThresh(),values[16]);
 assert.equal(configured.onlyHeavyAtomsForRms(),values[17]);
 assert.equal(configured.etVersion(),values[18]);
 assert.equal(configured.embedFragmentsSeparately(),values[19]);
 assert.equal(configured.useSmallRingTorsions(),values[20]);
 assert.equal(configured.useMacrocycleTorsions(),values[21]);
 assert.equal(configured.useMacrocycle14Config(),values[22]);
 assert.equal(configured.timeout(),values[23]);
 assert.equal(configured.cpci().size,0);
 assert.equal(configured.forceTransAmides(),values[25]);
 assert.equal(configured.useSymmetryForPruning(),values[26]);
 assert.equal(configured.boundsMatForceScaling(),values[27]);
 assert.equal(configured.trackFailures(),values[28]);
 assert.equal(configured.enableSequentialRandomSeeds(),values[29]);
 assert.equal(configured.symmetrizeConjugatedTerminalGroupsForPruning(),values[30]);

 for(const [index,values,type] of [[0,[-1,4294967296,1.5,NaN],RangeError],[1,[-2147483649,2147483648,Infinity],RangeError],[2,['42'],TypeError],[3,[1,'true'],TypeError]])for(const value of values){const args=Array(31).fill(undefined);args[index]=value;assert.throws(()=>new b.EmbedParams(...args),type);}const original=p.toJson(),changed=p.withJson('{"randomSeed":42,"trackFailures":true}');assert.equal(changed.randomSeed(),42);assert.equal(changed.trackFailures(),true);assert.equal(p.toJson(),original);assert.throws(()=>p.withJson('{'),e=>e.kind==='InvalidEmbedParametersJson'&&e.detail instanceof b.ConformerError&&typeof e.detail.detail==='string');assert.equal(p.toJson(),original);const getter=p.randomSeed;p.randomSeed=42;assert.equal(getter.call(p),-1);delete p.randomSeed;assert.equal(p.randomSeed(),-1);assert.equal(p.toJson(),original);
});
test('EmbedParams coordinate and CPCI maps retain null/empty distinction, signed keys and detached ordered values',()=>{
 const coordinates=new Map([[2,new Float64Array([4,5,6])],[-1,[1,2,3]]]),cpci=new Map([[[3,1],2.5],[[0,2],-1]]),args=Array(31).fill(undefined);args[8]=coordinates;args[24]=cpci;const p=new b.EmbedParams(...args);coordinates.get(-1)[0]=9;cpci.clear();const map=p.coordMap();assert.deepEqual([...map],[[-1,[1,2,3]],[2,[4,5,6]]]);map.get(-1)[0]=99;assert.equal(p.coordMap().get(-1)[0],1);assert.deepEqual([...p.cpci()],[[[0,2],-1],[[3,1],2.5]]);for(const bad of [{},new Map([[1,[1,2]]]),new Map([['1',[1,2,3]]]),new Map([[2147483648,[1,2,3]]])]){args[8]=bad;assert.throws(()=>new b.EmbedParams(...args));}args[8]=undefined;for(const bad of [new Map([[[1],2]]),new Map([[[0,-1],2]]),new Map([[[0,1],'2']])]){args[24]=bad;assert.throws(()=>new b.EmbedParams(...args));}
});
test('All sixteen embedding entrypoints preserve source-defined WASM SIGINT limitation and original receiver/parameters',()=>{
 const params=b.EmbedParams.dg().withJson('{"randomSeed":42,"trackFailures":true}'),original=params.toJson();for(const call of calls.filter(c=>c.params)){const m=b.Molecule.fromSmiles('C'),args=call.multiple?[2,params]:[params];assert.throws(()=>m[call.name](...args),e=>{assert.equal(e.kind,'Conformer');assert.ok(e.cause.detail instanceof b.ConformerRunError);assert.equal(e.cause.cause.message,'independent process SIGINT capability is unsupported on wasm32');assert.equal(e.message,'conformer operation failed: independent process SIGINT capability is unsupported on wasm32');return true;});assert.equal(m.num3dConformers(),0);assert.equal(params.toJson(),original);assert.equal(params.failures().length,0);}for(const [type,names] of [[b.EmbedMoleculeResult,['molecule','params','confId','ok']],[b.EmbedMultipleConfsResult,['molecule','params','confIds','requestedNumConfs','generatedCount']]])for(const name of names)assert.equal(typeof type.prototype[name],'function');
});
test('Default WASM failures keep original SIGINT source error, and supported invalid inputs remain failures with atomic receivers',()=>{
 for(const call of calls.filter(c=>!c.params)){const m=b.Molecule.fromSmiles('C');assert.throws(()=>m[call.name](...(call.multiple?[1]:[])),e=>{assert.equal(e.kind,'Conformer');assert.ok(e.cause.detail instanceof b.ConformerRunError);assert.equal(e.cause.cause.message,'independent process SIGINT capability is unsupported on wasm32');return true;});assert.equal(m.num3dConformers(),0);}
 const p=b.EmbedParams.dg().withJson('{"randomSeed":42}'),m=b.Molecule.fromSmiles('C');for(const n of [-1,4294967296,1.5,NaN])assert.throws(()=>m.with3dConformersWithParams(n,p),RangeError);assert.throws(()=>m.with3dConformersWithParams('1',p),TypeError);assert.equal(m.num3dConformers(),0);const bad=p.withJson('{"ETversion":3}');assert.throws(()=>m.embed3dConformerWithParams(bad),e=>e.kind==='Conformer'&&e.message.includes('Only version 1 and 2'));assert.equal(m.num3dConformers(),0);assert.equal(bad.etVersion(),3);const empty=b.Molecule.new();assert.throws(()=>empty.embed3dConformerWithParams(p),e=>e.kind==='Conformer'&&e.message.includes('molecule has no atoms'));assert.equal(empty.num3dConformers(),0);
});
test('Distance-geometry matrix query preserves nested numeric shape and leaves conformers unchanged',()=>{const m=b.Molecule.fromSmiles('CC'),matrix=m.dgBoundsMatrix();assert.equal(matrix.length,2);assert.equal(matrix[0].length,2);assert.equal(matrix[0][0],0);assert.equal(matrix[1][1],0);assert.ok(matrix[0][1]>matrix[1][0]);const old=matrix[0][1];matrix[0][1]=99;assert.equal(m.dgBoundsMatrix()[0][1],old);assert.equal(m.num3dConformers(),0);});
