import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });

test('SMILES grammar rejects repeated and dangling bond/separator tokens', () => {
    for (const smiles of ['C==C', 'C#', 'C..C', 'C\\\\\\C']) {
        assert.throws(() => ck.Molecule.fromSmiles(smiles), undefined, smiles);
    }
    const molecule = ck.Molecule.fromSmiles('C/C=C\\\\C');
    try {
        assert.equal(molecule.toSmiles(), 'C/C=C\\C');
    } finally {
        molecule.free();
    }
});

test('default MOL/SDF writing preserves stereo without source coordinates', () => {
    // Pinned RDKit 2026.03.1 canonical strings for the reported counterexamples.
    for (const [smiles, expected] of [
        ['N[C@@H](C)C(=O)O', 'C[C@H](N)C(=O)O'],
        ['N[C@H](C)C(=O)O', 'C[C@@H](N)C(=O)O'],
        ['F/C=C/F', 'F/C=C/F'],
        ['F/C=C\\F', 'F/C=C\\F'],
        ['C[C@H](O)[C@@H](O)C', 'C[C@H](O)[C@H](C)O'],
    ]) {
        const molecule = ck.Molecule.fromSmiles(smiles);
        try {
            const before = molecule.toBinary();
            for (const [write, read] of [['toMol', 'fromMol'], ['toSdf', 'fromSdf']]) {
                const restored = ck.Molecule[read](molecule[write]());
                try {
                    assert.equal(restored.toSmiles(), expected, `${smiles}: ${write}`);
                } finally {
                    restored.free();
                }
                assert.deepEqual(molecule.toBinary(), before, 'writer must preserve source');
            }
        } finally {
            molecule.free();
        }
    }
});

test('empty-molecule QED evaluates the source formula', () => {
    const molecule = ck.Molecule.fromSmiles('');
    try {
        assert.equal(molecule.qed(), 0.33942358984550913);
        assert.equal(molecule.toSmiles(), '');
    } finally {
        molecule.free();
    }
});
