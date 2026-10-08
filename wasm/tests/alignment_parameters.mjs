import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

const binding = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
binding.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });

test("all five alignment parameter constructors preserve source defaults and editable fields", () => {
    for (const name of ["AlignmentParameters", "BestAlignmentParameters", "CoordinateRmsdParameters", "AllConformerRmsdParameters", "ConformerAlignmentParameters"]) {
        const instance = new binding[name]();
        const factory = binding[name].new();
        for (const value of [instance, factory]) {
            assert.equal(value.weights, null);
            value.weights = [1, 2];
            const copy = value.weights;
            copy[0] = 9;
            assert.deepEqual(value.weights, [1, 2]);
            value.weights = null;
            assert.equal(value.weights, null);
            if ("probeConformerId" in value) {
                assert.equal(value.probeConformerId, -1);
                assert.equal(value.referenceConformerId, -1);
                value.probeConformerId = -2147483648;
                assert.equal(value.probeConformerId, -2147483648);
                assert.throws(() => { value.probeConformerId = 2147483648; }, RangeError);
            }
            if ("maxIterations" in value) {
                assert.equal(value.maxIterations, 50);
                assert.equal(value.reflect, false);
                value.maxIterations = 4294967295;
                assert.equal(value.maxIterations, 4294967295);
                assert.throws(() => { value.maxIterations = -1; }, RangeError);
            }
            if ("maxMatches" in value) {
                assert.equal(value.maxMatches, 1000000);
                assert.equal(value.symmetrizeConjugatedTerminalGroups, true);
                assert.deepEqual(value.atomMaps, []);
            }
            if ("ignoreHydrogens" in value) {
                assert.equal(value.ignoreHydrogens, true);
                assert.equal(value.numThreads, 1);
                value.numThreads = 0;
                assert.equal(value.numThreads, 0);
            }
            value.free();
        }
    }
});

test("atom maps and index arrays are detached without consuming input values", () => {
    const map = new binding.AlignmentAtomMap(3, 4);
    const params = binding.AlignmentParameters.new(-1, -1, [map], new Float64Array([2]), true, 99);
    assert.equal(map.probeAtom, 3);
    assert.equal(params.atomMap[0].referenceAtom, 4);
    map.probeAtom = 8;
    assert.equal(params.atomMap[0].probeAtom, 3);
    const detached = params.atomMap;
    detached[0].probeAtom = 11;
    assert.equal(params.atomMap[0].probeAtom, 3);
    assert.equal(params.reflect, true);
    assert.equal(params.maxIterations, 99);
    assert.deepEqual(params.weights, [2]);
    const best = new binding.BestAlignmentParameters(-1, -1, [[map]]);
    assert.equal(best.atomMaps[0][0].probeAtom, 8);
    assert.equal(map.referenceAtom, 4);
    const conformers = new binding.ConformerAlignmentParameters([0, 2], new Uint32Array([1, 4]));
    assert.deepEqual(conformers.atomIndices, [0, 2]);
    assert.deepEqual(conformers.conformerIds, [1, 4]);
    conformers.atomIndices = null;
    assert.equal(conformers.atomIndices, null);
    for (const n of [-1, 4294967296, 1.5, NaN, Infinity]) {
        assert.throws(() => binding.AlignmentAtomMap.new(n, 0), RangeError);
    }
    assert.throws(() => { params.weights = ["2"]; }, TypeError);
    assert.deepEqual(params.weights, [2]);
    assert.throws(() => { params.atomMap = [{ probeAtom: -1, referenceAtom: 0 }]; }, RangeError);
    assert.equal(params.atomMap[0].probeAtom, 3);
    for (const value of [map, params, best, conformers]) value.free();
});
