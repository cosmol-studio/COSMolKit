// Protocol samples only. Do not append a checklist of chemistry operations here.
const assert = require('node:assert/strict');
const ck = require('../pkg/cosmolkit_tsify_probe.js');

const p = new ck.RingSearchParams();
const m = new ck.Molecule('C1CC1');
assert.throws(() => { p.includeDativeBonds = 'false'; }, TypeError);
assert.equal(p.includeDativeBonds, false);
const configured = new ck.RingSearchParams({includeDativeBonds: true});
assert.equal(configured.includeDativeBonds, true);
assert.equal(configured.includeHydrogenBonds, p.includeHydrogenBonds);
assert.throws(() => new ck.RingSearchParams({includeDativeBonds: 'true'}), TypeError);
assert.throws(() => new ck.RingSearchParams({unknown: true}), TypeError);
const copied = configured.copy();
copied.includeDativeBonds = false;
assert.equal(configured.includeDativeBonds, true);
assert.throws(() => new ck.Molecule(6), TypeError);
assert.equal(m.symmetrizedSssr(p).numRings, 1);
assert.equal(m.symmetrizedSssr(p).numRings, 1); // borrowed Params is not consumed
assert.equal(m.symmetrizedSssr().numRings, 1);
assert.equal(m.symmetrizedSssr({includeDativeBonds: true}).numRings, 1);
assert.throws(() => m.symmetrizedSssr({unknown: true}), TypeError);
assert.throws(() => m.symmetrizedSssr(true), TypeError);
assert.throws(() => m.symmetrizedSssr(new ck.SmilesWriteParams()), TypeError);
const nested = new ck.NestedParams({child: p, optionalChild: p});
p.includeDativeBonds = true;
assert.equal(nested.child.includeDativeBonds, false); // construction by value
const child = nested.child;
child.includeHydrogenBonds = true;
assert.equal(nested.child.includeHydrogenBonds, true); // live child, not clone
const replacement = new ck.RingSearchParams({includeDativeBonds: true});
nested.child = replacement;
assert.equal(child.includeDativeBonds, true); // retained path resolves replacement
const optional = nested.optionalChild;
nested.optionalChild = undefined;
assert.throws(() => optional.includeDativeBonds, ReferenceError);
nested.free();
assert.equal(child.includeDativeBonds, true); // child owns state

// Shared native borrow: retained, readonly, and resolved after replacement.
const reaction = new ck.Reaction();
const view = reaction.matchParams();
assert.throws(() => { view.maxMatches = 17; }, TypeError);
const editable = view.copy();
editable.maxMatches = 17;
reaction.setMatchParams(editable);
assert.equal(view.maxMatches, 17);
reaction.free();
assert.equal(view.maxMatches, 17);

const selection = new ck.Selection({indices: [0, 1]});
const indices = selection.indices;
assert.throws(() => indices.push(2), TypeError);
selection.indices = [2];
assert.deepEqual(indices, [0, 1]);
assert.deepEqual(selection.indices, [2]);
assert.throws(() => { selection.indices = [-1]; });
assert.deepEqual(selection.indices, [2]); // failed write is atomic

// One map protocol witness: no writable receiver escapes through reads.
const layout = new ck.Coordinate2DParams({coordinateMap: new Map([[0, [1, 2]]])});
const map = layout.coordinateMap;
assert.deepEqual([...map], [[0, [1, 2]]]);
assert.equal(map.size, 1);
assert.throws(() => map.set(1, [2, 3]), TypeError);
assert.throws(() => Map.prototype.set.call(map, 1, [2, 3]), TypeError);
map.forEach((value, key, owner) => {
    assert.equal(owner, map);
    assert.throws(() => Map.prototype.clear.call(owner), TypeError);
});
assert.equal(map.valueOf(), map);
layout.coordinateMap = map; // read-only snapshots are valid field inputs
assert.deepEqual([...layout.coordinateMap], [[0, [1, 2]]]);

const input = {point: [1, 2, 3], indices: [0], mode: 'preserve'};
assert.deepEqual(ck.roundtrip(input), input);
assert.throws(() => ck.roundtrip({...input, point: [1, 2]}));
assert.throws(() => ck.roundtrip({...input, point: [1, 2, 3, 4]}), TypeError);
assert.throws(() => ck.roundtrip({...input, mode: 'unknown'}));
assert.throws(() => ck.roundtrip({...input, indices: [-1]}));
assert.deepEqual(ck.outcome({kind: 'completed', count: 2}), {kind: 'completed', count: 2});
assert.throws(() => ck.outcome({kind: 'completed', count: '2'}));
assert.deepEqual(ck.bytes(new Uint8Array([0, 255])), new Uint8Array([0, 255]));
assert.equal(ck.integer(9007199254740993n), 9007199254740993n);
// One numeric protocol witness, not a per-method validation inventory.
assert.equal(ck.Element.fromAtomicNumber(6).symbol, 'C');
assert.equal(ck.Element.fromAtomicNumber(255), undefined); // Valid u8, absent element.
for (const value of [262, -1, 6.9, NaN, Infinity, '6', true, null, undefined, {}]) {
    assert.throws(() => ck.Element.fromAtomicNumber(value), TypeError);
}
assert.throws(() => m.atom(2 ** 32), TypeError);
for (const value of [-1n, 1n << 64n, 6, '6']) {
    assert.throws(() => ck.integer(value), TypeError);
}
assert.throws(() => ck.bytes(new Uint16Array([262])), TypeError);

let saved;
const query = ck.QueryGraph.fromSmarts('C');
const accept = value => { saved = value; return true; };
const matcher = new ck.SubstructMatchParams({extraFinalCheck: accept});
assert.equal(matcher.extraFinalCheck, accept);
matcher.extraFinalCheck = undefined;
assert.equal(matcher.extraFinalCheck, undefined);
matcher.extraFinalCheck = accept;
matcher.extraFinalCheck = null;
assert.equal(matcher.extraFinalCheck, undefined);
matcher.extraFinalCheck = accept;
assert.equal(m.hasSubstructMatch(query, matcher), true);
const exception = new Error('user callback');
matcher.extraFinalCheck = () => { throw exception; };
assert.throws(() => m.hasSubstructMatch(query, matcher), e => e === exception);
matcher.extraFinalCheck = () => 1;
assert.throws(() => m.hasSubstructMatch(query, matcher), TypeError);
m.free();
assert.equal(saved.numAtoms, 3);
// Scoped native loans expire; only an explicit shared snapshot survives.
let loan, snapshot, candidate;
const tautomer = new ck.Molecule('CC(=O)C');
tautomer.enumerateTautomers({callback: (source, progress) => {
    loan = progress; snapshot = progress.snapshot(); candidate = progress.entry(0)[1];
    return true;
}});
assert.throws(() => loan.len, ReferenceError);
assert.throws(() => candidate.numAtoms, ReferenceError);
assert.ok(snapshot.entry(0)[1].numAtoms > 0);
assert.throws(() => tautomer.enumerateTautomers({callback: () => { throw exception; }}), e => e === exception);
assert.throws(() => ck.raise(), e => {
    assert.ok(e instanceof ck.OperationError && e instanceof Error);
    assert.equal(e.variant, ck.OperationErrorKind.Rings);
    const ring = e.cause;
    assert.equal(e.cause, ring);
    assert.ok(ring instanceof ck.RingFindingError && ring instanceof Error);
    assert.equal(ring.variant, ck.RingFindingErrorKind.InvalidTopology);
    const cause = ring.cause;
    assert.equal(ring.cause, cause);
    assert.ok(cause instanceof ck.TemplateAttachmentOrderError && cause instanceof Error);
    assert.equal(cause.variant, ck.TemplateAttachmentOrderErrorKind.TargetRemoved);
    assert.equal(cause.details.position, 2);
    assert.equal(cause.details.target, 3);
    assert.equal(e.details.rings.details.topology.details.atom, 1);
    return true;
});
// Native error projection shares the same Error protocol, without losing payload.
assert.throws(() => ck.Element.fromSymbol('not-an-element'), e => {
    assert.ok(e instanceof ck.ElementParseError && e instanceof Error);
    assert.equal(e.name, 'ElementParseError');
    assert.equal(e.message, "unknown element symbol 'not-an-element'");
    assert.equal(String(e), e.message);
    assert.equal(e.details.input, 'not-an-element');
    return true;
});
// One native iterator protocol witness; no per-collection inventory.
const iterator = ck.featureSpecs();
const sharedCursor = iterator[Symbol.iterator]();
const first = iterator.next();
assert.equal(first.done, false);
assert.ok(first.value instanceof ck.FeatureSpec);
iterator.free();
for (const value of sharedCursor) assert.ok(value instanceof ck.FeatureSpec);
assert.deepEqual(sharedCursor.next(), {done: true, value: undefined});
assert.deepEqual(sharedCursor.next(), {done: true, value: undefined});
console.log('runtime: 7 protocol groups passed');
