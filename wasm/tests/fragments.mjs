import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Fragment projections preserve ordering, last ties and the source', () => {
    const mol = ck.Molecule.fromSmiles('CC.OO.[Na+]');
    const before = mol.toSmiles();
    assert.deepEqual(mol.fragments().map(value => value.toSmiles()), ['CC', 'OO', '[Na+]']);
    assert.equal(mol.largestFragment().toSmiles(), 'OO');
    assert.equal(mol.toSmiles(), before);
    assert.deepEqual(ck.Molecule.new().fragments(), []);
    assert.throws(() => ck.Molecule.new().largestFragment(), error => error.kind === 'EmptyFragments');
});
