import assert from 'node:assert/strict';
import {test} from 'node:test';
import {pathToFileURL} from 'node:url';
import {readFileSync} from 'node:fs';
const b = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('MCS defaults, parameter/options forwarding and independent query results', () => {
    const inputs = ['CCO', 'CCN'].map(s => b.Molecule.fromSmiles(s));
    const before = inputs.map(m => m.toSmiles());
    const r = b.maximumCommonSubstructure(inputs);
    assert.deepEqual([r.atomCount, r.bondCount, r.completed, r.smarts], [2, 1, true, '[#6]-[#6]']);
    const query = r.query;
    assert.ok(query instanceof b.QueryGraph);
    assert.ok(inputs.every(m => m.hasSubstructMatch(query)));
    const p = new b.McsParameters();
    p.atomComparator = b.McsAtomComparator.Any;
    for (const opts of [p, {atomComparator: b.McsAtomComparator.Any}]) {
        const result = b.maximumCommonSubstructure(inputs, opts);
        assert.deepEqual([result.atomCount, result.bondCount, result.completed], [3, 2, true]);
        assert.equal(result.smarts, '[#6]-[#6]-[#8,#7]');
    }
    assert.equal(b.maximumCommonSubstructureWithParams(inputs, p).atomCount, 3);
    assert.deepEqual(inputs.map(m => m.toSmiles()), before);
    inputs.forEach(m => m.free());
    assert.equal(query.numAtoms(), 2);
});

test('MCS nested configuration views are live and failed assignment is atomic', () => {
    const inputs = ['[NH4+]', 'N'].map(s => b.Molecule.fromSmiles(s));
    assert.equal(b.maximumCommonSubstructure(inputs).atomCount, 1);
    const p = new b.McsParameters();
    p.atomCompareParameters.matchFormalCharge = true;
    p.bondCompareParameters.matchStereo = true;
    p.timeout = 10;
    assert.equal(p.atomCompareParameters.matchFormalCharge, true);
    assert.equal(p.bondCompareParameters.matchStereo, true);
    assert.equal(b.maximumCommonSubstructure(inputs, p).atomCount, 0);
    assert.equal(b.maximumCommonSubstructure(inputs, {atomCompareParameters: {matchFormalCharge: true}}).atomCount, 0);
    assert.throws(() => { p.timeout = -1; }, RangeError);
    assert.equal(p.timeout, 10);
    assert.throws(() => b.maximumCommonSubstructure(inputs, {unknownOption: true}), TypeError);
    assert.throws(() => b.maximumCommonSubstructure(inputs, {atomCompareParameters: {unknownOption: true}}), TypeError);
});

test('MCS empty, store-all and typed error paths', () => {
    const empty = b.maximumCommonSubstructure(['Cl', 'Br'].map(s => b.Molecule.fromSmiles(s)));
    assert.deepEqual([empty.atomCount, empty.bondCount, empty.query, empty.smarts, empty.completed], [0, 0, null, '', true]);
    const inputs = ['CCO', 'CCN'].map(s => b.Molecule.fromSmiles(s));
    const stored = b.maximumCommonSubstructure(inputs, {storeAll: true});
    assert.ok(stored.degenerate instanceof Map);
    assert.ok(stored.degenerate.size > 0);
    assert.throws(() => b.maximumCommonSubstructure([]), e => e.name === 'McsError' && e.domain === 'mcs' && e.kind === 'State' && e.detail instanceof b.McsError);
    assert.throws(() => b.maximumCommonSubstructure(inputs, {threshold: 1.1}), e => e.name === 'McsError');
    assert.throws(() => b.maximumCommonSubstructure([{}]), TypeError);
});
