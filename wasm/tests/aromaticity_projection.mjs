import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";
const b = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });

test("aromaticity enum, frozen parameters and all four operations project canonical semantics", () => {
    assert.equal(new b.AromaticityParams().model, b.AromaticityModel.Rdkit);
    for (const model of [b.AromaticityModel.Rdkit, b.AromaticityModel.Simple, b.AromaticityModel.Mdl, b.AromaticityModel.Mmff94]) {
        const p = new b.AromaticityParams(model);
        assert.equal(p.model, model);
        assert.throws(() => { p.model = b.AromaticityModel.Custom; }, TypeError);
        const source = b.Molecule.fromSmilesWithSanitize("C1=CC=CC=C1", false);
        const before = source.toSmiles();
        const value = source.withAssignedAromaticityWithParams(p);
        // Pinned Aromaticity.cpp setMMFFAromaticity requires SP2 carbon/nitrogen.
        // This fixture explicitly disables sanitization/hybridization assignment.
        assert.equal(value.toSmiles(), model === b.AromaticityModel.Mmff94 ? "C1=CC=CC=C1" : "c1ccccc1");
        assert.equal(source.toSmiles(), before);
        assert.equal(source.assignAromaticityWithParams(p), undefined);
        assert.equal(source.toSmiles(), value.toSmiles());
        value.free(); source.free(); p.free();
    }
    const source = b.Molecule.fromSmilesWithSanitize("C1=CC=CC=C1", false);
    const before = source.toSmiles(), value = source.withAssignedAromaticity();
    assert.equal(value.toSmiles(), "c1ccccc1");
    assert.equal(source.toSmiles(), before);
    assert.equal(source.assignAromaticity(), undefined);
    assert.equal(source.toSmiles(), value.toSmiles());
    value.free(); source.free();
    assert.throws(() => new b.AromaticityParams(77));
});

test("source-defined unsupported custom model keeps structured cause and failed state", () => {
    const source = b.Molecule.fromSmilesWithSanitize("C1=CC=CC=C1", false);
    const before = source.toSmiles(), params = new b.AromaticityParams(b.AromaticityModel.Custom);
    for (const invoke of [() => source.withAssignedAromaticityWithParams(params), () => source.assignAromaticityWithParams(params)]) {
        assert.throws(invoke, error => {
            assert.ok(error instanceof Error);
            assert.equal(error.domain, "operation");
            assert.equal(error.kind, "Aromaticity");
            assert.equal(error.cause.name, "AromaticityError");
            assert.equal(error.cause.domain, "aromaticity");
            assert.equal(error.cause.kind, "UnsupportedModel");
            assert.equal(error.cause.model, b.AromaticityModel.Custom);
            assert.ok(error.cause.reason.length > 0);
            assert.ok(error.cause.detail instanceof b.AromaticityError);
            assert.equal(error.cause.detail.kind, "UnsupportedModel");
            assert.equal(error.cause.detail.cause, null);
            return true;
        });
        assert.equal(source.toSmiles(), before);
    }
    params.free(); source.free();
});
