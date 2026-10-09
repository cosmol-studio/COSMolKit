import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });
const preset = process.env.COSMOLKIT_WASM_PRESET;

test("preset exports exactly its selected optional domains", () => {
    const molecule = ck.Molecule.fromSmiles("CCO");
    try {
        assert.equal(molecule.numAtoms(), 3);
        assert.equal(molecule.toSmiles(), "CCO");
        const search = ["full", "core-search", "core-analysis", "core-reaction"].includes(preset);
        assert.equal(typeof molecule.substructMatches, search ? "function" : "undefined");
        if (search) {
        const query = ck.QueryGraph.fromSmarts("[#6]");
        try {
            assert.equal(molecule.substructMatches(query).length, 2);
        } finally {
            query.free();
        }
        }
        for (const [method, enabled] of [
            ["fingerprintMorgan", ["full", "core-fingerprints", "core-analysis"].includes(preset)],
            ["qed", ["full", "core-analysis"].includes(preset)],
            ["with3dConformer", ["full", "core-3d"].includes(preset)],
            ["toBinary", false],
            ["toSvg", ["full", "core-depict"].includes(preset)],
            ["toInchi", ["full", "core-inchi"].includes(preset)],
        ]) {
            assert.equal(typeof molecule[method], enabled ? "function" : "undefined", `${preset}: ${method}`);
        }
        assert.equal(typeof ck.Reaction, ["full", "core-reaction"].includes(preset) ? "function" : "undefined");
        assert.equal(typeof ck.BioStructure, ["full", "core-bio"].includes(preset) ? "function" : "undefined");
        assert.equal(typeof ck.Molecule.fromBinary, "undefined");
        assert.equal(typeof ck.PickleError, "undefined");
        const batch = ck.MoleculeBatch.fromSmilesList(["CCO"]);
        try {
            assert.deepEqual(batch.toSmilesList(), ["CCO"]);
            assert.equal(typeof batch.toSvgList, ["full", "core-depict"].includes(preset) ? "function" : "undefined");
            assert.equal(typeof batch.dgBoundsMatrixList, ["full", "core-3d"].includes(preset) ? "function" : "undefined");
        } finally { batch.free(); }
    } finally {
        molecule.free();
    }
});

test("ordinary MOL/SDF IO works without search or depiction, without capability leakage", () => {
    const molecule = ck.Molecule.fromXyzBlock("2\nethane fragment\nC 0 0 0\nC 1.5 0 0\n");
    const params = new ck.MolBlockWriteParams(undefined, undefined, false);
    try {
        const text = molecule.toSdfWithParams(params);
        const record = ck.SdfRecord.fromSdf(text);
        const restored = record.molecule();
        try {
            assert.equal(restored.numAtoms(), 2);
            const search = ["full", "core-search", "core-analysis", "core-reaction"].includes(preset);
            // The detached value is always registered in Rust as metadata.
            // Only its parser and search operations are optional capabilities.
            assert.equal(typeof ck.QueryGraph, "function");
            assert.equal(typeof ck.QueryGraph.fromSmarts, search ? "function" : "undefined");
            assert.equal(typeof ck.parseSmarts, search ? "function" : "undefined");
            assert.equal(typeof record.queryGraph, "function");
            assert.throws(() => record.queryGraph(), /graph|query/i);
        } finally { restored.free(); record.free(); }
    } finally { params.free(); molecule.free(); }
});

test("missing depiction reports an explicit capability error, never silent stereo loss", () => {
    const molecule = ck.Molecule.fromSmiles("C[C@H](F)Cl");
    try {
        const before = molecule.toSmiles();
        if (["full", "core-depict"].includes(preset)) {
            const restored = ck.Molecule.fromMol(molecule.toMol());
            try { assert.equal(restored.toSmiles(), before); }
            finally { restored.free(); }
        } else {
            assert.throws(() => molecule.toMol(), error =>
                error.kind === "MolWrite" && error.cause.kind === "MissingCapability"
                && error.cause.capability === "depict");
            assert.throws(() => molecule.toSdf2d(), /depict/);
        }
        assert.equal(molecule.toSmiles(), before);
        assert.equal(molecule.coordinates2d().length, 0);
    } finally { molecule.free(); }
});
