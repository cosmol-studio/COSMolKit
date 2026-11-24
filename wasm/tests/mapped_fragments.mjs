import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as ck from '../pkg/index.js';

test('component molecules and source indices survive their WASM handles', () => {
    const source = ck.Molecule.fromSmiles('C1.CC.C1');
    const [first, second] = source.fragmentsWithAtomMapping({ sanitize: false });
    assert.deepEqual([first.sourceAtomIndices, second.sourceAtomIndices], [[0, 3], [1, 2]]);

    const owned = first.molecule;
    const cut = ck.BondCut.atBond(0);
    owned.cutBonds_([cut]);
    assert.equal(owned.numAtoms, 4);
    assert.equal(source.numAtoms, 4);
    source.cutBonds_([cut]);
    source.free();
    owned.free();
    cut.free();

    const indices = first.sourceAtomIndices;
    const molecule = first.molecule;
    first.free();
    second.free();
    assert.deepEqual(indices, [0, 3]);
    assert.ok(Object.isFrozen(indices));
    assert.equal(molecule.toSmiles(), 'CC');
    molecule.free();
});

test('fragment sanitization uses the requested setting', () => {
    const source = ck.Molecule.fromSmiles('C.C(F)(F)(F)(F)F', { sanitize: false, removeHs: false });
    assert.throws(() => source.fragmentsWithAtomMapping(), error =>
        error instanceof ck.OperationError && error.cause instanceof ck.MoleculeFragmentsError &&
        error.cause.details.componentIndex === 1);

    const config = new ck.FragmentExtractionParams({ sanitize: false });
    for (const parts of [source.fragmentsWithAtomMapping({ sanitize: false }),
                         source.fragmentsWithAtomMapping(config), source.fragmentsWithAtomMappingWithParams(config)]) {
        assert.deepEqual(parts.map(p => p.sourceAtomIndices), [[0], [1, 2, 3, 4, 5, 6]]);
        const molecule = parts[1].molecule;
        assert.equal(molecule.numBonds, 5);
        molecule.free();
        for (const part of parts) part.free();
    }
    config.free();
    source.free();
});
