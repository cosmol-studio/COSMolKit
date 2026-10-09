import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('Avalon writable configuration and typed failures retain causes', () => {
    const mol = ck.Molecule.fromSmiles('CCO'), before = mol.toSmiles();
    const params = new ck.AvalonFingerprintParams();
    params.nBits = 7; params.isQuery = false; params.bitFlags = 32767;
    assert.match(params.toString(), /nBits=7/);
    for (const config of [params, {nBits: 7}]) {
        assert.throws(() => mol.fingerprintAvalon(config), error => error.kind === 'InvalidArguments'
            && error.cause.name === 'AvalonEngineError'
            && error.cause.detail instanceof ck.AvalonEngineError);
    }
    assert.throws(() => { params.isQuery = 1; }, TypeError);
    assert.equal(mol.toSmiles(), before);
});
test('Avalon conversion succeeds with coordinates without requiring depict', () => {
    const block = '\n  RDKit          2D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n'
        + '   -1.2990   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
        + '    0.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
        + '    1.2990   -0.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n'
        + '  1  2  1  0\n  2  3  1  0\nM  END\n';
    const mol = ck.Molecule.fromMol(block);
    assert.deepEqual(Array.from(mol.fingerprintAvalon({nBits:64}).onBits()), [6,14,30,31,42]);
});
