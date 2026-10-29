import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('dense fingerprint values preserve checked inputs, detached bits and typed errors',()=>{
 const input=new Uint32Array([7,1,1]),a=b.Fingerprint.fromOnBits(8,input);input[0]=0;
 assert.equal(a.nBits(),8);assert.deepEqual(a.onBits(),[1,7]);const bits=a.onBits();bits.push(0);assert.deepEqual(a.onBits(),[1,7]);
 const other=b.Fingerprint.fromOnBits(8,[1,2]);assert.equal(a.tanimoto(other),1/3);assert.equal(a.tanimoto(a),1);assert.equal(other.nBits(),8);
 assert.deepEqual(b.Fingerprint.fromOnBits(0,[]).onBits(),[]);
 assert.throws(()=>b.Fingerprint.fromOnBits(8,[8]),e=>{assert.equal(e.domain,'fingerprints');assert.equal(e.kind,'SparseIndexOutOfRange');assert.equal(e.index,8n);assert.equal(e.size,8n);assert.ok(e.detail instanceof b.FingerprintError);assert.equal(e.detail.cause,null);return true;});
 assert.throws(()=>a.tanimoto(b.Fingerprint.fromOnBits(7,[])),e=>e.kind==='BitLengthMismatch'&&e.left===8n&&e.right===7n);
 for(const n of [-1,1.5,2**32,NaN,Infinity]){assert.throws(()=>b.Fingerprint.fromOnBits(n,[]),RangeError);assert.throws(()=>b.Fingerprint.fromOnBits(8,[n]),RangeError);}
 assert.throws(()=>b.Fingerprint.fromOnBits('8',[]),TypeError);assert.throws(()=>b.Fingerprint.fromOnBits(8,{}),TypeError);
});
for(const [Class,wide] of [[b.SparseCountFingerprint,true],[b.SparseCountFingerprint32,false]])test(`${Class.name} transports exact widths, detached maps and all canonical arithmetic`,()=>{
 const key=n=>wide?BigInt(n):n,a=Class.new(key(8)),other=Class.new(key(8));assert.equal(a.length(),key(8));assert.equal(a.value(key(0)),0);
 a.setValue(key(1),2);a.setValue(key(3),4);other.setValue(key(1),3);other.setValue(key(2),5);
 assert.deepEqual([...a.nonzeroElements()],[[key(1),2],[key(3),4]]);const map=a.nonzeroElements();map.set(key(1),99);assert.equal(a.value(key(1)),2);
 assert.equal(a.totalValue(),6);assert.equal(a.totalValue(true),6);assert.equal(a.fuzzyAnd(other).value(key(1)),2);assert.equal(a.fuzzyOr(other).value(key(1)),3);
 assert.equal(a.withAdded(other).value(key(1)),5);assert.equal(a.withSubtracted(other).value(key(1)),-1);
 assert.equal(a.withAddedScalar(2).value(key(1)),4);assert.equal(a.withSubtractedScalar(2).value(key(1)),0);assert.equal(a.withMultipliedScalar(3).value(key(1)),6);assert.equal(a.withDividedScalar(2).value(key(3)),2);
 assert.equal(a.withAdded(a).value(key(1)),4);assert.equal(a.value(key(1)),2);assert.equal(other.value(key(1)),3);
 a.setValue(key(3),0);assert.equal(a.nonzeroElements().has(key(3)),false);a.setValue(key(1),-2);assert.equal(a.totalValue(),-2);assert.equal(a.totalValue(true),2);
 assert.throws(()=>a.value(key(8)),e=>e.kind==='SparseIndexOutOfRange'&&e.index===8n&&e.size===8n);
 assert.throws(()=>a.fuzzyAnd(Class.new(key(9))),e=>e.kind==='BitLengthMismatch');
 assert.throws(()=>a.withDividedScalar(0),e=>e.kind==='UndefinedArithmetic'&&typeof e.site==='string');
 for(const n of [-2147483649,2147483648,NaN,1.5])assert.throws(()=>a.setValue(key(0),n),RangeError);
 assert.throws(()=>a.totalValue('false'),TypeError);
 const maximum=wide?(1n<<64n)-1n:2**32-1,huge=Class.new(maximum);huge.setValue(maximum,7);assert.equal(huge.value(maximum),7);assert.equal(huge.length(),maximum);
 if(wide){const big=2n**53n+1n;huge.setValue(big,-3);assert.equal(huge.nonzeroElements().get(big),-3);for(const n of [-1n,1n<<64n])assert.throws(()=>Class.new(n),RangeError);for(const v of [8,'8',null])assert.throws(()=>Class.new(v),TypeError);}
 else for(const n of [-1,2**32,1.5,NaN])assert.throws(()=>Class.new(n),RangeError);
});
