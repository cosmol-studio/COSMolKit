import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });

test("concrete fingerprints and configured providers do not require public search", () => {
    const molecule = ck.Molecule.fromSmiles("CCO");
    const generator = ck.MorganFingerprintGenerator.new();
    const featureProvider = ck.MorganAtomInvariantsGenerator.features([]);
    const configured = new ck.MorganFingerprintGenerator(new ck.MorganParams(0), featureProvider);
    try {
        assert.deepEqual(molecule.fingerprintMorgan().onBits(), molecule.fingerprintMorganWithGenerator(generator).onBits());
        assert.ok(molecule.fingerprintMorganWithGenerator(configured) instanceof ck.Fingerprint);
        const params = new ck.TopologicalFingerprintParams();
        const request = new ck.TopologicalFingerprintOutputRequest(true, true);
        const output = molecule.fingerprintTopologicalWithOutputWithParams(params, request);
        assert.deepEqual(output.fingerprint().onBits(), molecule.fingerprintTopological().onBits());
        assert.equal(output.atomBits().length, 3);
        assert.ok(output.bitInfo().size > 0);
        const mask = ck.Fingerprint.fromOnBits(2048, [674]);
        const seed = new ck.LayeredFingerprintParams(undefined, undefined, undefined, undefined, [10, 20, 30], mask);
        const layered = molecule.fingerprintLayeredWithOutputWithParams(seed);
        assert.deepEqual(layered.fingerprint().onBits(), [674]);
        assert.deepEqual(layered.atomCounts(), [11, 22, 31]);
        assert.deepEqual(molecule.fingerprintPattern().onBits(), molecule.fingerprintPatternWithParams(new ck.PatternFingerprintParams()).onBits());
        assert.throws(() => molecule.fingerprintPatternWithParams(new ck.PatternFingerprintParams(0)), error => error.kind === "EmptyFingerprint");
        assert.equal(molecule.toSmiles(), "CCO");
    } finally { configured.free(); featureProvider.free(); generator.free(); molecule.free(); }
});
