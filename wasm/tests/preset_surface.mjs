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
        const query = ck.QueryGraph.fromSmarts("[#6]");
        try {
            assert.equal(molecule.substructMatches(query).length, 2);
        } finally {
            query.free();
        }
        for (const [method, enabled] of [
            ["morganFingerprint", ["full", "core-fingerprints", "core-analysis"].includes(preset)],
            ["qed", ["full", "core-analysis"].includes(preset)],
            ["with3dConformer", ["full", "core-3d"].includes(preset)],
            ["toBinary", preset === "full"],
            ["toInchi", ["full", "core-inchi"].includes(preset)],
        ]) {
            assert.equal(typeof molecule[method], enabled ? "function" : "undefined", `${preset}: ${method}`);
        }
        assert.equal(typeof ck.Reaction, ["full", "core-reaction"].includes(preset) ? "function" : "undefined");
        assert.equal(typeof ck.BioStructure, ["full", "core-bio"].includes(preset) ? "function" : "undefined");
    } finally {
        molecule.free();
    }
});
