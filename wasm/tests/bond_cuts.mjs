import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as ck from '../pkg/index.js';

test('cut methods preserve attachment sides and differ only in ownership', () => {
    const source = new ck.Molecule('CCO');
    const cuts = [ck.BondCut.atBond(1, 101, 202)];
    const product = source.withBondsCut(cuts);
    assert.equal(source.toSmiles(), 'CCO');
    assert.deepEqual(product.fragments().map(f => f.toSmiles()), ['[101*]CC', '[202*]O']);
    assert.equal(source.cutBonds_(cuts), undefined);
    assert.equal(source.toSmiles(), product.toSmiles());
});

test('late label overflow preserves the molecule and identifies the failing cut', () => {
    const source = new ck.Molecule('CCO');
    const cuts = [ck.BondCut.atBond(0), ck.BondCut.atBond(1, 65536, 8)];
    for (const method of ['withBondsCut', 'cutBonds_']) {
        assert.throws(() => source[method](cuts), e => {
            assert.ok(e instanceof ck.OperationError);
            assert.ok(e.cause instanceof ck.BondCutError);
            const error = e.details.bondCut;
            assert.equal(error.variant, ck.BondCutErrorKind.IsotopeOutOfRange);
            assert.deepEqual([error.details.bond, error.details.attachment, error.details.label], [1, 1, 65536]);
            return true;
        });
        assert.equal(source.toSmiles(), 'CCO');
    }
    assert.throws(() => ck.BondCut.atBond(0, undefined, 0), e =>
        e instanceof ck.BondCutError && e.variant === ck.BondCutErrorKind.IncompleteLabels);
});
