import assert from "node:assert/strict";
import test from "node:test";
import {readFileSync} from "node:fs";
import {pathToFileURL} from "node:url";
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const identifier = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3";
const key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N";

test("InChI defaults, configurable calls, typed errors and reusable parameters", () => {
    const molecule = ck.Molecule.fromSmiles("CCO");
    const read = new ck.InchiReadParams();
    const write = new ck.InchiWriteParams();
    assert.equal(molecule.toInchi(), identifier);
    assert.equal(molecule.toInchiKey(), key);
    assert.equal(ck.inchiToInchiKey(identifier), key);
    write.options = "-AuxNone";
    for (const params of [write, {options: "-AuxNone"}]) {
        assert.equal(molecule.toInchi(params), identifier);
        assert.equal(molecule.toInchiKey(params), key);
    }
    assert.equal(molecule.toInchiWithParams(write), identifier);
    for (const sanitize of [false, true]) for (const remove of [false, true]) {
        read.sanitize = sanitize;
        read.removeHydrogens = remove;
        for (const params of [read, {sanitize, removeHydrogens: remove}]) {
            const parsed = ck.Molecule.fromInchi(identifier, params);
            assert.equal(parsed.toInchi(), identifier);
            parsed.free();
        }
        const parsed = ck.Molecule.fromInchiWithParams(identifier, read);
        assert.equal(parsed.toInchi(), identifier);
        parsed.free();
    }
    for (const call of [() => ck.Molecule.fromInchi("not InChI"), () => ck.inchiToInchiKey("not InChI")]) {
        assert.throws(call, e => e instanceof Error && e.name === "InchiError"
            && typeof e.kind === "number" && e.domain === "inchi" && typeof e.operation === "string"
            && e.detail instanceof ck.InchiError);
    }
    assert.throws(() => molecule.toInchi({bogus: true}), TypeError);
    assert.throws(() => molecule.toInchi({options: 3}), TypeError);
    assert.throws(() => ck.Molecule.fromInchi(identifier, {sanitize: "false"}), TypeError);
    assert.throws(() => { read.sanitize = "true"; }, TypeError);
    assert.throws(() => { write.options = 3; }, TypeError);
    assert.equal(write.options, "-AuxNone");
    assert.equal(molecule.toSmiles(), "CCO");
    write.free(); read.free(); molecule.free();
});
