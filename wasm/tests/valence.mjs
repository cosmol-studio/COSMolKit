import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

const b = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });
const raw = text => b.Molecule.fromSmilesWithParams(text,
    new b.SmilesParseParams(false, undefined, undefined, undefined, false));
const missing = e => e.name === "AtomPairReadError" && e.kind === "Preparation"
    && e.cause.kind === "MissingPreparedValence";
const invalidValence = e => e instanceof Error && e.name === "OperationError"
    && e.kind === "Valence" && e.cause.name === "ValenceError"
    && e.cause.kind === "InvalidValence" && e.cause.atom === 0
    && e.cause.atomicNumber === 6 && e.cause.detail instanceof b.ValenceError;

test("Valence parameters preserve canonical defaults and reject host coercion", () => {
    const params = new b.ValenceParams();
    assert.equal(params.model, b.ValenceModel.RdkitLike);
    assert.equal(params.strict, true);
    assert.equal(new b.ValenceParams(b.ValenceModel.RdkitLike, false).strict, false);
    for (const value of [-1, 1, 1.5, NaN, Infinity]) {
        assert.throws(() => new b.ValenceParams(value), RangeError);
    }
    for (const value of [null, "0", false, {}]) {
        assert.throws(() => new b.ValenceParams(value), TypeError);
    }
    for (const value of [null, "false", 0, {}]) {
        assert.throws(() => new b.ValenceParams(undefined, value), TypeError);
    }
    assert.throws(() => { params.strict = false; }, TypeError);
});

test("All four valence operations retain value semantics and typed failure atomicity", () => {
    for (const explicit of [false, true]) {
        const value = raw("CCO"), params = new b.ValenceParams();
        assert.throws(() => value.atomPairFingerprint(), missing);
        const prepared = explicit ? value.withAssignedValenceWithParams(params)
            : value.withAssignedValence();
        assert.ok(prepared.atomPairFingerprint() instanceof b.Fingerprint);
        assert.throws(() => value.atomPairFingerprint(), missing);
        assert.equal(explicit ? value.assignValenceWithParams(params) : value.assignValence(), undefined);
        assert.ok(value.atomPairFingerprint() instanceof b.Fingerprint);
        value.free();
        assert.equal(prepared.toSmiles(), "CCO");

        const invalid = raw("C(F)(F)(F)(F)F"), before = invalid.toSmiles();
        assert.throws(() => explicit ? invalid.withAssignedValenceWithParams(params)
            : invalid.withAssignedValence(), invalidValence);
        assert.throws(() => explicit ? invalid.assignValenceWithParams(params)
            : invalid.assignValence(), invalidValence);
        assert.equal(invalid.toSmiles(), before);
        assert.throws(() => invalid.atomPairFingerprint(), missing);
        const relaxed = new b.ValenceParams(undefined, false);
        assert.ok(invalid.withAssignedValenceWithParams(relaxed).atomPairFingerprint() instanceof b.Fingerprint);
        assert.throws(() => invalid.atomPairFingerprint(), missing);
        invalid.assignValenceWithParams(relaxed);
        assert.ok(invalid.atomPairFingerprint() instanceof b.Fingerprint);
    }
    for (const params of [null, undefined, {}, new b.SmilesParseParams()]) {
        assert.throws(() => raw("CC").withAssignedValenceWithParams(params));
        assert.throws(() => raw("CC").assignValenceWithParams(params));
    }
});

test("The sanitize convenience parser retains structured SMILES errors", () => {
    for (const sanitize of [true, false]) {
        assert.throws(() => b.Molecule.fromSmilesWithSanitize("[", sanitize), e =>
            e instanceof Error && e.name === "SmilesError" && e.detail instanceof b.SmilesError
            && e.kind === "Parse" && e.cause instanceof Error);
    }
    for (const sanitize of [undefined, null, "false", 0, 1, {}]) {
        assert.throws(() => b.Molecule.fromSmilesWithSanitize("CCO", sanitize), TypeError);
    }
});
