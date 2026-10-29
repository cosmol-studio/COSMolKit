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

const constraints = [
    ['Distance', [0, 1], [false, 1, 1.1, 20], {relative: false, minDistance: 1, maxDistance: 1.1, forceConstant: 20}],
    ['Angle', [0, 1, 2], [false, 80, 90, 20], {relative: false, minAngleDeg: 80, maxAngleDeg: 90, forceConstant: 20}],
    ['Torsion', [0, 1, 2, 3], [false, -10, 10, 20], {relative: false, minDihedralDeg: -10, maxDihedralDeg: 10, forceConstant: 20}],
    ['Position', [0], [0.01, 20], {maxDisplacement: 0.01, forceConstant: 20}],
];
for (const model of ['UFF', 'MMFF94', 'MMFF94s']) {
    for (const [kind, atoms, values, options] of constraints) {
        test(`${model} ${kind} constraint instances, required options, setters and failure atomicity`, () => {
            const source = ck.Molecule.fromSmiles('CCCC').withHydrogens().with3dConformer({randomSeed: 42});
            const factory = () => model === 'UFF' ? source.uffForceField() : source.mmffForceField({mmffVariant: model});
            const field = factory();
            const peer = factory();
            const untouched = factory();
            const original = Array.from(source.coordinates3d(0));
            const Class = ck[`ForceField${kind}ConstraintParams`];
            const params = new Class(...values);
            const method = `add${kind}Constraint_`;
            assert.throws(() => new Class(), TypeError);
            for (const key of Object.keys(options)) {
                const missing = {...options};
                delete missing[key];
                assert.throws(() => peer[method](...atoms, missing), TypeError);
            }
            assert.throws(() => peer[method](...atoms, {...options, unknown: 1}), TypeError);
            assert.throws(() => { params.forceConstant = 'wrong'; }, TypeError);
            assert.equal(params.forceConstant, 20);
            params.forceConstant = 30;
            field[method](...atoms, params);
            peer[method](...atoms, {...options, forceConstant: 30});
            const displaced = field.position(0);
            displaced[0] += 0.2;
            field.setPosition_(0, displaced);
            peer.setPosition_(0, displaced);
            assert.equal(field.energy(), peer.energy());
            assert.deepEqual(field.gradient(), peer.gradient());
            const before = field.energy();
            const invalid = [...atoms];
            invalid[0] = 999;
            assert.throws(() => field[method](...invalid, params), error => {
                assert.equal(error.kind, ck.MolecularForceFieldErrorKind.InvalidAtomIndex);
                assert.equal(error.atomIndex, 999);
                assert.ok(error.detail instanceof ck.ForceFieldError);
                return true;
            });
            assert.equal(field.energy(), before);
            assert.deepEqual(Array.from(source.coordinates3d(0)), original);
            assert.deepEqual(untouched.position(0), Array.from(original.slice(0, 3)));
            for (const value of [params, field, peer, untouched, source]) value.free();
        });
    }
}
