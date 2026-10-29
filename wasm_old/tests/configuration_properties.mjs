import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import test from 'node:test';

const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('parser assignment uses constructor validation and preserves other fields', () => {
    const params = new ck.SmilesParseParams();
    params.parseName = false;
    params.sanitize = false;
    assert.equal(params.parseName, false);
    const configured = ck.Molecule.fromSmiles('C(F)(F)(F)(F)F', {sanitize: false});
    assert.equal(configured.numAtoms(), 6);
    configured.free();
    assert.throws(() => ck.Molecule.fromSmiles('CCO', {sanitise: false}), TypeError);
    assert.equal(params.sanitize, false);
    assert.throws(() => { params.sanitize = 'false'; }, TypeError);
    assert.equal(params.sanitize, false);
    assert.equal(params.parseName, false);
    const molecule = ck.Molecule.fromSmilesWithParams('C(F)(F)(F)(F)F', params);
    assert.equal(molecule.numAtoms(), 6);
    assert.throws(() => ck.Molecule.fromSmiles('C(F)(F)(F)(F)F'));
    molecule.free();
    params.free();
});

test('hydrogen assignment is used by the actual operation and isolates arrays', () => {
    const params = new ck.AddHsParams();
    const selection = [2];
    params.onlyOnAtoms = selection;
    params.addCoords = false;
    selection[0] = 0;
    assert.deepEqual(params.onlyOnAtoms, [2]);
    assert.throws(() => { params.onlyOnAtoms = [-1]; }, RangeError);
    assert.deepEqual(params.onlyOnAtoms, [2]);
    const molecule = ck.Molecule.fromSmiles('CCO');
    const result = molecule.withHydrogensWithParams(params);
    assert.equal(result.numAtoms(), 4);
    for (const configuration of [params, {onlyOnAtoms: [2]}]) {
        const same = molecule.withHydrogens(configuration);
        assert.equal(same.numAtoms(), result.numAtoms());
        assert.equal(same.toSmiles(), result.toSmiles());
        same.free();
    }
    assert.throws(() => molecule.withHydrogens({onlyOnAtoms: [-1]}), RangeError);
    assert.throws(() => molecule.withHydrogens({unknown: true}), TypeError);
    assert.equal(molecule.numAtoms(), 3);
    result.free();
    molecule.free();
    params.free();
});

test('batch configuration validates enum assignment atomically', () => {
    const params = new ck.BatchParams();
    params.nJobs = 1;
    params.errors = ck.BatchErrorMode.KeepErrors;
    assert.equal(params.nJobs, 1);
    assert.equal(params.errors, ck.BatchErrorMode.KeepErrors);
    assert.throws(() => { params.errors = 99; }, RangeError);
    assert.equal(params.errors, ck.BatchErrorMode.KeepErrors);
    assert.equal(params.nJobs, 1);
    params.errors = null;
    assert.equal(params.errors, null);
    params.free();
});
