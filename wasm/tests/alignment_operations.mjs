import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";
const binding = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
binding.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });

function source(shift) {
    const row = (x, symbol) => `${x.toFixed(4).padStart(10)}    0.0000    0.0000 ${symbol}   0  0  0  0  0  0  0  0  0  0  0  0`;
    return binding.Molecule.fromSdf(`alignment\n     RDKit          3D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n${row(shift, "C")}\n${row(shift + 1, "O")}\n  1  2  1  0\nM  END\n$$$$\n`);
}

test("all eighteen alignment operations expose canonical results and mutation semantics", () => {
    const probe = source(3), reference = source(0);
    const params = new binding.AlignmentParameters();
    const best = new binding.BestAlignmentParameters();
    const coordinate = new binding.CoordinateRmsdParameters();
    const all = new binding.AllConformerRmsdParameters();
    const conformers = new binding.ConformerAlignmentParameters();
    const before = [...probe.coordinates3d(0)];
    for (const result of [probe.alignmentTransformTo(reference), probe.alignmentTransformToWithParams(reference, params), probe.bestAlignmentTo(reference), probe.bestAlignmentToWithParams(reference, best)]) {
        assert.ok(result instanceof binding.AlignmentResult);
        assert.ok(result.rmsd() < 1e-12);
        const transform = result.transform();
        const matrix = transform.matrix();
        assert.deepEqual(matrix.map(row => row.length), [4, 4, 4, 4]);
        assert.ok(Math.abs(matrix[0][3] + 3) < 1e-12);
        matrix[0][3] = 999;
        assert.ok(Math.abs(transform.matrix()[0][3] + 3) < 1e-12);
        const map = result.atomMap();
        assert.equal(map.length, 2);
        assert.equal(map[0].probeAtom, 0);
        map[0].probeAtom = 55;
        assert.equal(result.atomMap()[0].probeAtom, 0);
    }
    assert.ok(probe.bestRmsdTo(reference) < 1e-12);
    assert.ok(probe.bestRmsdToWithParams(reference, best) < 1e-12);
    assert.equal(probe.coordinateRmsdTo(reference), 3);
    assert.equal(probe.coordinateRmsdToWithParams(reference, coordinate), 3);
    assert.deepEqual(probe.allConformerBestRmsds(), []);
    assert.deepEqual(probe.allConformerBestRmsdsWithParams(all), []);
    for (const [value, result] of [probe.withAlignmentTo(reference), probe.withAlignmentToWithParams(reference, params)]) {
        assert.ok(result.rmsd() < 1e-12);
        assert.ok(value.coordinateRmsdTo(reference) < 1e-12);
        assert.deepEqual([...probe.coordinates3d(0)], before);
    }
    for (const [value, report] of [probe.withAlignedConformers(), probe.withAlignedConformersWithParams(conformers)]) {
        assert.ok(value instanceof binding.Molecule);
        assert.ok(report instanceof binding.ConformerAlignmentReport);
        assert.deepEqual([...report.rmsds()], []);
    }
    assert.ok(probe.alignTo(probe).rmsd() < 1e-12);
    assert.ok(probe.alignToWithParams(probe, params).rmsd() < 1e-12);
    assert.ok(probe.alignTo(reference).rmsd() < 1e-12);
    assert.ok(probe.alignToWithParams(reference, params).rmsd() < 1e-12);
    assert.ok(probe.coordinateRmsdTo(reference) < 1e-12);
    assert.deepEqual([...probe.alignConformers().rmsds()], []);
    assert.deepEqual([...probe.alignConformersWithParams(conformers).rmsds()], []);
    assert.equal(reference.numAtoms(), 2);
});

test("alignment errors preserve domain, kind, fields, cause and wrapper empty-list policy", () => {
    const probe = source(3), reference = source(0);
    const params = new binding.AlignmentParameters(-1, -1, [], []);
    assert.ok(probe.alignmentTransformToWithParams(reference, params).rmsd() < 1e-12);
    assert.deepEqual(params.atomMap, []);
    assert.deepEqual(params.weights, []);
    params.weights = [1];
    const before = [...probe.coordinates3d(0)];
    assert.throws(() => probe.alignmentTransformToWithParams(reference, params), error => {
        assert.ok(error instanceof Error);
        assert.equal(error.name, "AlignmentError");
        assert.equal(error.domain, "alignment");
        assert.equal(error.kind, "WeightCountMismatch");
        assert.equal(error.mapLen, 2);
        assert.equal(error.weightLen, 1);
        assert.ok(error.detail instanceof binding.AlignmentError);
        assert.equal(error.detail.kind, error.kind);
        return true;
    });
    assert.throws(() => probe.alignToWithParams(reference, params), error => {
        assert.equal(error.domain, "operation");
        assert.equal(error.kind, "Alignment");
        assert.equal(error.cause.domain, "alignment");
        assert.equal(error.cause.kind, "WeightCountMismatch");
        assert.equal(error.cause.mapLen, 2);
        return true;
    });
    assert.deepEqual([...probe.coordinates3d(0)], before);
    params.weights = null;
    params.probeConformerId = 77;
    assert.throws(() => probe.alignmentTransformToWithParams(reference, params), error => {
        assert.equal(error.kind, "ConformerNotFound");
        assert.equal(error.id, 77);
        return true;
    });
    assert.throws(() => binding.Molecule.fromSmiles("CO").alignmentTransformTo(reference), error => error.kind === "NoConformers");
});
