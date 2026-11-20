import assert from 'node:assert/strict';
import { test } from 'node:test';
import * as ck from '../pkg/index.js';

test('BRICS sites preserve rule orientation and repeated fragments', () => {
    const mol = new ck.Molecule('CCCOCCC');
    const sites = mol.findBricsBonds();
    assert.deepEqual(sites.map(s => [s.bond, s.ends]), [[2, [[3, 3], [2, 4]]], [3, [[3, 3], [4, 4]]]]);
    assert.equal(mol.withBondsCut([sites[0].toBondCut()]).toSmiles(), '[3*]OCCC.[4*]CCC');
    assert.equal(mol.cutBonds_(sites.map(s => s.toBondCut())), undefined);
    assert.deepEqual(mol.fragments().map(f => f.toSmiles()), ['[4*]CCC', '[3*]O[3*]', '[4*]CCC']);
});

test('regular r-BRICS labels are usable with the existing cut operation', () => {
    const mol = new ck.Molecule('CCCCCCCC');
    const [site] = mol.findRBricsBonds();
    assert.ok(site instanceof ck.RBricsBondMatch);
    assert.deepEqual(site.ends, [[3, 20], [4, 21]]);
    assert.deepEqual(mol.decomposeRBrics().sort(), ['[20*]CCCC', '[21*]CCCC']);
    mol.cutBonds_([site.toBondCut()]);
    assert.equal(mol.toSmiles(), '[20*]CCCC.[21*]CCCC');
});

test('recursive BRICS returns unique terminal SMILES or the product error', () => {
    const mol = new ck.Molecule('CCCOCCC');
    assert.deepEqual(mol.decomposeBrics().sort(), ['[3*]O[3*]', '[4*]CCC']);
    assert.equal(mol.toSmiles(), 'CCCOCCC');
    const complex = new ck.Molecule('O=C([O-])C1=CC=C(N2N=C3[N]4=C2C5=CC=CC=C5[O-][Fe+3]4[O-]C6=CC=CC=C63)C=C1');
    assert.throws(() => complex.decomposeBrics(), e =>
        e instanceof ck.BricsDecompositionError && e.variant === ck.BricsDecompositionErrorKind.UnkeyedProduct &&
        e.cause instanceof ck.SanitizeError);
});

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

test('RECAP retains shared products and distinguishes missing keys from leaves', () => {
    const g = new ck.Molecule('CC(=O)Nc1ccc(NC(=O)C)cc1').decomposeRecap();
    const intermediate = '*Nc1ccc(NC(C)=O)cc1';
    assert.deepEqual(g.parents('*C(C)=O').sort(), [intermediate, g.root]);
    assert.deepEqual(g.children(g.root).sort(), ['*C(C)=O', intermediate]);
    assert.deepEqual(g.leaves().sort(), ['*C(C)=O', '*Nc1ccc(N*)cc1']);
    assert.deepEqual(g.children('*C(C)=O'), []);
    for (const method of ['parents', 'children']) {
        assert.throws(() => g[method](42), TypeError);
        assert.throws(() => g[method]('absent'), e =>
            e instanceof ck.RecapLookupError && e.details.product === 'absent');
    }
    assert.deepEqual(new ck.Molecule('C').decomposeRecap().leaves(), []);
});
