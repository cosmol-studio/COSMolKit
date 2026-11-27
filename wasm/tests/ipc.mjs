import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as packaged from '../pkg/index.js';
import * as raw from '../pkg/cosmolkit_wasm.js';

for (const [name, api] of [['packaged', packaged], ['raw', raw]]) {
    test(`${name}: Ipc and AvgIpc match the exact RDKit build through the last bit`, () => {
        const bits = value => new BigUint64Array(Float64Array.of(value).buffer)[0];
        // b9518c8a: the five-cycle isolates entropy rounding; the larger graph
        // accumulates it over several nonzero coefficients.
        for (const [smiles, ipc, avgIpc] of [
            ['C1CCCC1', 0x4036e2ed6df50345n, 0x3ffc2afcd619dca4n],
            ['O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3', 0x410f098c5f48ecb7n, 0x40070c1a44280aedn],
        ]) {
            const mol = api.Molecule.fromSmiles(smiles);
            try {
                assert.equal(bits(mol.ipc()), ipc);
                assert.equal(bits(mol.avgIpc()), avgIpc);
            } finally { mol.free(); }
        }
    });

    test(`${name}: matrix overrides preserve values and matrix ownership`, () => {
        const mol = api.Molecule.fromSmiles('CC');
        const disconnected = api.Molecule.fromSmiles('C.C');
        const matrix = disconnected.distanceMatrix();
        try {
            assert.equal(mol.ipc(), 2);
            assert.equal(mol.avgIpc(), 1);
            assert.equal(mol.ipc(matrix), 0);
            assert.equal(mol.avgIpcWithDistanceMatrix(matrix), 0);
            assert.equal(mol.ipcWithDistanceMatrix([[0, 1], [0, 0]]), 0);
            assert.equal(mol.avgIpc([[0, NaN], [1, 0]]), 1);
            assert.equal(mol.avgIpc(matrix), 0);
            assert.throws(() => mol.ipc([[0, 1], [0]]), error => {
                assert.equal(error.variant, api.IpcErrorKind.InvalidDistanceMatrixShape);
                const details = error.details;
                try { assert.deepEqual([details.row, details.actual, details.expected], [1, 1, 2]); }
                finally { details.free(); error.free(); }
                return true;
            });
            assert.throws(() => mol.ipc([[0, 0], ['1', 0]]), TypeError);
            assert.equal(mol.ipc(), 2);
        } finally {
            matrix.free(); disconnected.free(); mol.free();
        }
    });
}

test('out-of-range and overflow return JavaScript numbers', () => {
    const params = { sanitize: false, removeHs: false };
    const path = packaged.Molecule.fromSmilesWithParams('C'.repeat(1473), params);
    const outOfRange = raw.Molecule.fromSmilesWithParams('C'.repeat(4000), params);
    try {
        assert.equal(path.ipc(), Infinity);
        assert.ok(Number.isNaN(outOfRange.avgIpc()));
    } finally {
        path.free(); outOfRange.free();
    }
});
