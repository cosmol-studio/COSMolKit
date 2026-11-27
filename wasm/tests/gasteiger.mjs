import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as packaged from '../pkg/index.js';
import * as raw from '../pkg/cosmolkit_wasm.js';

for (const [name, api] of [['packaged', packaged], ['raw', raw]]) {
    test(`${name}: charges preserve graph rows, and invalid iterations never coerce`, () => {
        const mol = api.Molecule.fromSmiles('[Fe]C.CCO');
        const zero = new api.GasteigerChargeParams(0);
        let rows = [];
        try {
            for (const value of [-1, 1.25, NaN, Infinity, 2 ** 31]) {
                assert.throws(() => new api.GasteigerChargeParams(value));
                assert.throws(() => { zero.iterations = value; });
                assert.throws(() => mol.gasteigerCharges({ iterations: value }));
                assert.equal(zero.iterations, 0);
            }

            for (const invalid of [{iterations: null}, {throwOnParamFailure: null}, {unknown: true}]) {
                assert.throws(() => mol.gasteigerCharges(invalid), TypeError);
            }
            for (const configured of [zero, {iterations: 0, throwOnParamFailure: undefined}, zero]) {
                rows = mol.gasteigerCharges(configured);
                assert.ok(rows.every(row => row.atomCharge === 0 && row.attachedHydrogenCharge === 0));
                rows.splice(0).forEach(row => row.free());
            }
            rows = mol.gasteigerChargesWithParams(zero);
            assert.ok(rows.every(row => row.atomCharge === 0 && row.attachedHydrogenCharge === 0));
            rows.splice(0).forEach(row => row.free());
            rows = mol.gasteigerCharges({iterations: undefined});
            assert.equal(rows[4].atomCharge, -0.3966637147812598);
            rows.splice(0).forEach(row => row.free());
            rows = mol.gasteigerCharges();
            assert.ok(Number.isNaN(rows[0].atomCharge));
            assert.equal(rows[0].attachedHydrogenCharge, 0);
            assert.ok(Number.isNaN(rows[1].attachedHydrogenCharge));
            assert.equal(rows[4].atomCharge, -0.3966637147812598);
            assert.equal(rows[4].attachedHydrogenCharge, 0.21002230656786322);
            assert.ok(Number.isNaN(mol.minPartialCharge()));
            assert.equal(mol.peoeVsa14(), mol.peoeVsa()[13]);
            assert.throws(() => mol.gasteigerCharges({ throwOnParamFailure: true }));
        } finally {
            rows.forEach(row => row.free());
            zero.free();
            mol.free();
        }
    });
}
