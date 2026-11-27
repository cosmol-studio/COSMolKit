import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as api from '../pkg/index.js';

test('BCUT2D returns all eight ethanol extrema or a typed bond failure', () => {
    const ethanol = api.Molecule.fromSmiles('CCO');
    const dative = api.Molecule.fromSmiles('[2H]C([2H])N->B');
    const reference = [16.248991329327914, 10.908168035086371, 1.3402935186707783,
        -1.52235869276019, 1.2962426382141046, -1.5683635470355741,
        3.8358023355894861, 0.31840561317391147];
    try {
        const values = ethanol.bcut2d();
        assert.equal(values.length, reference.length);
        reference.forEach((r, i) => assert.ok(Math.abs(values[i] - r) <= 1e-12 * Math.max(1, Math.abs(r))));
        assert.throws(() => dative.bcut2d(), err => {
            const cause = err.cause;
            const details = cause.details;
            try {
                assert.equal(cause.variant, api.DescriptorErrorKind.UnsupportedBcutBond);
                assert.equal(details.bond, 3);
                assert.equal(details.order, api.BondOrder.Dative);
                return true;
            } finally {
                details.free(); cause.free(); err.free();
            }
        });
    } finally {
        ethanol.free(); dative.free();
    }
});
