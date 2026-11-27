import assert from 'node:assert/strict';
import { test } from 'node:test';
import { Molecule } from '../pkg/index.js';

test('EState vectors, projection direction and graph edits', () => {
    const mol = Molecule.fromSmiles('CCO');
    try {
        const held = mol.estateIndices();
        const area = mol.estateVsa();
        const indices = mol.vsaEstate();
        assert.equal(area.length, 11);
        assert.equal(indices.length, 10);
        assert.ok(Math.abs(area[5] - 6.923737199690624) < 1e-12);
        assert.equal(mol.estateVsa6(), area[5]);
        assert.equal(mol.vsaEstate3(), held[2]);
        assert.equal(indices[2], held[2]);

        mol.removeAtom_(1);
        assert.deepEqual([...mol.estateIndices()], [0, 0]);
        assert.equal(mol.estateVsa6(), 0);
        assert.ok(Math.abs(held[2] - 545 / 72) < 1e-12);
    } finally {
        mol.free();
    }
});
