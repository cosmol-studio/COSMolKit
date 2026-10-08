import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

test("Element and elementInfo preserve canonical value and result boundaries", async () => {
    const binding = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
    binding.initSync({ module: readFileSync(process.env.COSMOLKIT_WASM_BINARY) });
    const carbon = binding.Element.fromSymbol("C");
    assert.equal(carbon.atomicNumber(), 6);
    assert.equal(carbon.symbol(), "C");
    assert.equal(carbon.atomicNumber(), 6);
    assert.equal(binding.Element.fromSymbol("Uut").symbol(), "Nh");
    assert.equal(binding.Element.fromSymbol("Uup").symbol(), "Mc");
    for (const invalid of ["", "c", " C", "C\0", "é"]) {
        assert.equal(binding.Element.fromSymbol(invalid), null);
    }
    assert.equal(binding.Element.fromAtomicNumber(119), null);
    assert.equal(binding.Element.fromAtomicNumber(255), null);
    for (const invalid of [-1, 256, 6.5, NaN, Infinity, -Infinity]) {
        assert.throws(() => binding.Element.fromAtomicNumber(invalid), RangeError);
    }
    const info = binding.elementInfo(carbon);
    assert.equal(info.element().atomicNumber(), 6);
    assert.equal(info.symbol(), "C");
    assert.equal(info.atomicNumber(), 6);
    assert.equal(info.period(), 2);
    assert.equal(info.outerElectrons(), 4);
    assert.deepEqual([...info.valences()], [4]);
    assert.equal(info.rb0(), 0.77);
    assert.equal(info.atomicWeight(), 12.011);
    const values = info.valences();
    values[0] = -123;
    assert.deepEqual([...info.valences()], [4]);
    assert.equal(carbon.atomicNumber(), 6);
    info.free();
    carbon.free();
});
