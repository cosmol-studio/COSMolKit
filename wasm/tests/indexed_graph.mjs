import assert from 'node:assert/strict';
import {test} from 'node:test';
import * as ck from '../pkg/index.js';

test('canonical order, bond selection and coordinates use the same source rows', () => {
    const source = ck.Molecule.fromSmiles('OCC');
    source.add3dConformer_([[1, 0, 0], [2, 0, 0], [3, 0, 0]]);
    const ranks = source.canonicalAtomRanks();
    const order = ranks.map((_, i) => i).sort((i, j) => ranks[i] - ranks[j]);
    assert.deepEqual(order, [2, 0, 1]);
    const ordered = source.withAtomsRenumbered(order);
    const path = ordered.shortestPath(0, 2);
    assert.deepEqual(path, [0, 2]);
    const bond = ordered.bondBetweenAtoms(...path);
    assert.equal(bond, 1);
    assert.equal(ordered.bondBetweenAtoms(0, 1), undefined);
    const matrix = ordered.adjacencyMatrix();
    assert.equal(matrix.dimension, 3);
    assert.equal(matrix.get(0, 2), 1);
    matrix.free();
    const selection = ck.SubgraphSelection.bonds([bond]);
    const result = ordered.extractSubgraph(selection);
    const retained = ordered.retainSubgraph_(selection);
    const owned = result.molecule, other = result.molecule, mapping = result.mapping;
    assert.deepEqual(retained.sourceAtomIndices, [0, 2]);
    assert.deepEqual(retained.sourceBondIndices, [1]);
    assert.deepEqual(owned.coordinates3d().map(row => row[0]), [3, 2]);
    assert.deepEqual(source.coordinates3d().map(row => row[0]), [1, 2, 3]);
    for (const handle of [source, selection, retained, result, ordered]) handle.free();
    owned.removeAtom_(0);
    owned.free();
    assert.equal(other.numAtoms, 2);
    assert.deepEqual(mapping.sourceBondIndices, [1]);
    mapping.free(); other.free();
});

test('failed in-place edits preserve the receiver; insertion uses an independent handle', () => {
    const m = ck.Molecule.fromSmiles('OCC');
    const env = m.atomEnvironment(1, 0);
    assert.deepEqual(env.atomDistances, new Map([[0, 0], [1, 1]]));
    env.free();
    assert.deepEqual(m.subgraphsInRange(1, 2)[1], [[0, 1]]);
    const invalid = ck.SubgraphSelection.atomsAndBonds([0], [0]);
    assert.throws(() => m.retainSubgraph_(invalid), e => e.cause instanceof ck.IndexedGraphError);
    invalid.free();
    assert.throws(() => m.renumberAtoms_([0, 1, 0]), ck.OperationError);
    assert.equal(m.numAtoms, 3);
    const clone = m.withMoleculeAppended(m);
    assert.throws(() => m.insertMolecule_(m), TypeError);
    assert.equal(m.numAtoms, 3);
    const right = m.copy();
    try { m.insertMolecule_(right); }
    finally { right.free(); }
    assert.equal(m.toSmiles(), clone.toSmiles());
    assert.equal(m.bondBetweenAtoms(3, 4), 2);
    assert.deepEqual(m.shortestPath(0, 3), []);
    clone.free(); m.free();
});
