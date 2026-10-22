// Protocol samples only. Do not append a checklist of chemistry operations here.
const assert = require('node:assert/strict');
const ck = require('../pkg/cosmolkit_tsify_probe.js');

const p = new ck.RingSearchParams();
const m = new ck.Molecule('C1CC1');
assert.equal(m.num_rings(p), 1);
assert.equal(m.num_rings(p), 1); // borrowed Params is not consumed
const nested = new ck.NestedParams(p);
p.include_dative_bonds = true;
assert.equal(nested.child.include_dative_bonds, false); // construction by value
const child = nested.child;
child.include_hydrogen_bonds = true;
assert.equal(nested.child.include_hydrogen_bonds, true); // live child, not clone
nested.free();
assert.equal(child.include_hydrogen_bonds, true); // child owns state

const selection = new ck.Selection([0, 1]);
const indices = selection.indices;
assert.throws(() => indices.push(2), TypeError);
selection.indices = [2];
assert.deepEqual(indices, [0, 1]);
assert.deepEqual(selection.indices, [2]);
assert.throws(() => { selection.indices = [-1]; });
assert.deepEqual(selection.indices, [2]); // failed write is atomic

const input = {point: [1, 2, 3], indices: [0], mode: 'preserve'};
assert.deepEqual(ck.roundtrip(input), input);
assert.throws(() => ck.roundtrip({...input, point: [1, 2]}));
assert.throws(() => ck.roundtrip({...input, mode: 'unknown'}));
assert.throws(() => ck.roundtrip({...input, indices: [-1]}));
assert.deepEqual(ck.outcome({kind: 'completed', count: 2}), {kind: 'completed', count: 2});
assert.throws(() => ck.outcome({kind: 'completed', count: '2'}));
assert.deepEqual(ck.bytes(new Uint8Array([0, 255])), new Uint8Array([0, 255]));
assert.equal(ck.integer(9007199254740993n), 9007199254740993n);

let saved;
assert.equal(m.inspect(value => { saved = value; return true; }), true);
const exception = new Error('user callback');
assert.throws(() => m.inspect(() => { throw exception; }), e => e === exception);
assert.throws(() => m.inspect(() => 1), TypeError);
m.free();
assert.equal(saved.num_atoms(), 3);
assert.throws(() => ck.raise({kind: 'Syntax', offset: 3, message: 'input'}, exception), e => {
    assert.ok(e instanceof Error);
    assert.equal(e.cause, exception);
    assert.deepEqual(ck.failure_details(e), {kind: 'Syntax', offset: 3, message: 'input'});
    return true;
});
console.log('runtime: 6 protocol groups passed');
