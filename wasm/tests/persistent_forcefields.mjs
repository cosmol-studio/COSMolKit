import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import test from 'node:test';
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

for (const flavor of ['mmff', 'uff']) {
    test(`${flavor} persistent evaluator owns coordinates, pins, derivatives and two-step minimization`, () => {
        const source = ck.Molecule.fromSmiles('CCO').withHydrogens().with3dConformer({randomSeed: 42});
        const original = Array.from(source.coordinates3d(0));
        const field = source[`${flavor}ForceField`]();
        const Class = flavor === 'mmff' ? ck.MmffForceFieldParams : ck.UffForceFieldParams;
        const other = source[`${flavor}ForceField`](new Class());
        const configured = source[`${flavor}ForceField`]({conformerId: 0});
        assert.equal(configured.energy(), field.energy());
        configured.free();
        const reference = source[`${flavor}EnergyGradient`]();
        assert.equal(field.energy(), reference.energy());
        assert.deepEqual(field.gradient().flat(), Array.from(reference.gradient()));
        reference.free();
        const positions = field.positions();
        positions[1][0] += 0.2;
        for (const handle of [field, other]) {
            handle.setPositions_(positions);
            handle.setFixedAtoms_([0]);
            assert.deepEqual(handle.fixedAtoms(), [0]);
            assert.deepEqual(handle.gradient()[0], [0, 0, 0]);
            assert.ok(handle.gradientUnconstrained().flat().every(Number.isFinite));
            const result = handle.energyGradient();
            assert.equal(result.energy(), handle.energy());
            assert.deepEqual(result.gradient(), handle.gradient());
            result.free();
        }
        assert.deepEqual(Array.from(source.coordinates3d(0)), original);
        source.free();
        const pinned = field.position(0);
        const before = field.energy();
        const options = new ck.ForceFieldMinimizeParams(2);
        options.energyTolerance = 1e-6;
        assert.deepEqual([options.maxIterations, options.forceTolerance, options.energyTolerance], [2, 1e-4, 1e-6]);
        assert.throws(() => { options.maxIterations = -1; }, RangeError);
        assert.equal(options.maxIterations, 2);
        const outcome = field.minimize_(options);
        const explicit = other.minimizeWithParams_(options);
        assert.ok(outcome instanceof ck.ForceFieldMinimizeOutcome);
        assert.ok(outcome.iterations() <= 2);
        assert.ok(outcome.energy() <= before);
        assert.equal(outcome.energy(), explicit.energy());
        assert.deepEqual(field.positions(), other.positions());
        assert.deepEqual(field.position(0), pinned);
        const stable = field.positions();
        assert.throws(() => field.setPositions_([]), error => {
            assert.equal(error.kind, ck.MolecularForceFieldErrorKind.CoordinateCount);
            assert.ok(error.detail instanceof ck.ForceFieldError);
            assert.equal(error.detail.kind(), ck.MolecularForceFieldErrorKind.CoordinateCount);
            assert.equal(error.detail.actual(), 0);
            assert.equal(error.detail.expected(), stable.length);
            return true;
        });
        assert.throws(() => field.setPosition_(0, [NaN, 0, 0]), error => error.kind === ck.MolecularForceFieldErrorKind.NonFiniteCoordinate && error.atomIndex === 0 && error.component === 0);
        assert.throws(() => field.setFixedAtoms_([999]), error => error.kind === ck.MolecularForceFieldErrorKind.InvalidFixedAtom);
        assert.deepEqual(field.positions(), stable);
        assert.deepEqual(field.fixedAtoms(), [0]);
        for (const value of [field, other, outcome, explicit, options]) value.free();
    });
}

test('persistent factories expose typed missing-conformer errors and checked option fields', () => {
    const molecule = ck.Molecule.fromSmiles('CCO');
    for (const flavor of ['mmff', 'uff']) {
        assert.throws(() => molecule[`${flavor}ForceField`](), error => {
            assert.equal(error.kind, ck.MolecularForceFieldErrorKind.MissingConformer);
            assert.ok(error.detail instanceof (flavor === 'mmff' ? ck.MmffForceFieldError : ck.UffForceFieldError));
            assert.equal(error.detail.requested(), null);
            return true;
        });
        assert.throws(() => molecule[`${flavor}ForceField`]({conformerId: -1}), RangeError);
        assert.throws(() => molecule[`${flavor}ForceField`]({unknown: 1}), TypeError);
    }
    molecule.free();
});
