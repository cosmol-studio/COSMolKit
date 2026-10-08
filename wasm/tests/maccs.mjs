import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('All MACCS calls retain source raw and public projections and detached values',()=>{
 const p=new b.MaccsFingerprintParams();assert.equal(p.nBits,166);assert.throws(()=>{p.nBits=64;},TypeError);
 for(const text of ['CCO','c1ccncc1','CC(=O)N','CC.O']){const m=b.Molecule.fromSmiles(text),before=m.toSmiles(),fp=m.maccsFingerprint(),raw=m.maccsFingerprintRaw();assert.ok(fp instanceof b.Fingerprint);assert.equal(fp.nBits(),166);assert.equal(raw.nBits(),167);assert.deepEqual(fp.onBits(),m.maccsFingerprintWithParams(p).onBits());assert.deepEqual(fp.onBits(),raw.onBits().map(v=>v-1));assert.equal(m.toSmiles(),before);m.free();assert.equal(fp.nBits(),166);assert.ok(fp.onBits().length>0);}
});
test('Unsupported MACCS widths preserve exact source category, option and reason',()=>{
 const m=b.Molecule.fromSmiles('CCO');for(const width of [0,64,167]){const p=new b.MaccsFingerprintParams(width);assert.equal(p.nBits,width);assert.throws(()=>m.maccsFingerprintWithParams(p),e=>e.name==='MaccsFingerprintError'&&e.kind==='UnsupportedOption'&&e.domain==='Fingerprint'&&e.option==='MaccsFingerprintParams.n_bits'&&typeof e.reason==='string'&&e.detail instanceof b.MaccsFingerprintError&&e.detail.option===e.option&&e.detail.reason===e.reason&&e.detail.bit===null&&e.detail.cause===null);}
 assert.throws(()=>new b.MaccsFingerprintParams(-1),RangeError);assert.throws(()=>new b.MaccsFingerprintParams(1.5),RangeError);assert.throws(()=>new b.MaccsFingerprintParams('166'),TypeError);assert.equal(m.toSmiles(),'CCO');
});
