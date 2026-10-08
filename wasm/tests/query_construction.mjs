import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('complete query construction routes preserve detached graph counts and names',()=>{
 const p=new b.SmartsParseParams();assert.deepEqual([p.allowCxsmiles,p.strictCxsmiles,p.parseName,p.mergeHs,p.skipCleanup,p.debugParse],[true,true,true,false,false,false]);assert.deepEqual(p.replacements,new Map());
 for(const q of [b.QueryGraph.fromSmarts('[#6]-[#8] alcohol'),b.QueryGraph.fromSmartsWithParams('[#6]-[#8] alcohol',p),b.parseSmarts('[#6]-[#8] alcohol'),b.parseSmartsWithParams('[#6]-[#8] alcohol',p)]){assert.ok(q instanceof b.QueryGraph);assert.equal(q.numAtoms(),2);assert.equal(q.numBonds(),1);assert.equal(q.name(),'alcohol');}
 const replacements=new Map([['{X}','[#8]']]),configured=new b.SmartsParseParams(false,false,false,true,true,true,replacements);replacements.set('{X}','[#7]');assert.deepEqual(configured.replacements,new Map([['{X}','[#8]']]));assert.deepEqual([configured.allowCxsmiles,configured.strictCxsmiles,configured.parseName,configured.mergeHs,configured.skipCleanup,configured.debugParse],[false,false,false,true,true,true]);
 const simple=new b.SmartsParseParams(undefined,undefined,undefined,undefined,undefined,undefined,{'{X}':'[#8]'});assert.equal(b.QueryGraph.fromSmartsWithParams('{X}',simple).numAtoms(),1);assert.equal(b.QueryGraph.fromSmarts('C').name(),null);const copy=simple.replacements;copy.clear();assert.equal(simple.replacements.size,1);
 assert.throws(()=>{p.mergeHs=true;},TypeError);assert.throws(()=>new b.SmartsParseParams('true'),TypeError);assert.throws(()=>new b.SmartsParseParams(undefined,undefined,undefined,undefined,undefined,undefined,{X:3}),TypeError);
});
test('invalid SMARTS preserves typed parser category and source position',()=>{
 for(const call of [()=>b.QueryGraph.fromSmarts('['),()=>b.QueryGraph.fromSmartsWithParams('[',new b.SmartsParseParams()),()=>b.parseSmarts('['),()=>b.parseSmartsWithParams('[',new b.SmartsParseParams())])assert.throws(call,e=>e instanceof Error&&e.name==='SmartsParseError'&&e.domain==='search'&&e.kind==='UnclosedBracket'&&e.position===0&&e.detail instanceof b.SmartsParseError);
 assert.throws(()=>b.QueryGraph.fromSmarts('C1'),e=>e.kind==='Parse'&&e.message==='SMARTS parse error: unclosed ring');
});
