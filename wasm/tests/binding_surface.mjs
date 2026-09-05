import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

test("WASM binding preserves runtime values and errors", async () => {
    const modulePath = process.env.COSMOLKIT_WASM_MODULE;
    const wasmPath = process.env.COSMOLKIT_WASM_BINARY;
    assert.ok(modulePath, "COSMOLKIT_WASM_MODULE must point to generated wasm-bindgen JS");
    assert.ok(wasmPath, "COSMOLKIT_WASM_BINARY must point to generated wasm-bindgen binary");

    const binding = await import(pathToFileURL(modulePath).href);
    binding.initSync({ module: readFileSync(wasmPath) });

    const molecule = binding.Molecule.fromSmiles("CCO");
    assert.equal(molecule.toSmiles(), "CCO");
    assert.equal(molecule.numAtoms(), 3);
    assert.equal(molecule.numBonds(), 2);
    assert.equal(molecule.molecularWeight(), 46.069);
    assert.ok(Math.abs(molecule.exactMolecularWeight() - 46.041864812) < 1e-9);
    assert.ok(Number.isFinite(molecule.crippenLogP()));
    assert.ok(Number.isFinite(molecule.crippenMolarRefractivity()));
    assert.ok(Math.abs(molecule.tpsa() - 20.23) < 1e-9);
    assert.equal(molecule.hBondAcceptors(), 1);
    assert.equal(molecule.hBondDonors(), 1);
    assert.equal(molecule.formula(), "C2H6O");
    assert.equal(molecule.numHeavyAtoms(), 3);
    assert.equal(molecule.numRings(), 0);
    assert.equal(molecule.numAromaticRings(), 0);
    assert.equal(molecule.numRotatableBonds(), 0);
    assert.ok(Number.isFinite(molecule.qed()));
    assert.ok(Number.isFinite(molecule.hallKierAlpha()));
    assert.ok(Number.isFinite(molecule.kappa1()));
    assert.ok(Number.isFinite(molecule.kappa2()));
    assert.ok(Number.isFinite(molecule.kappa3()));
    assert.ok(Number.isFinite(molecule.chi0()));
    assert.ok(Number.isFinite(molecule.chi1()));
    assert.ok(Number.isFinite(molecule.chi(1, true, false)));
    assert.ok(Number.isFinite(molecule.phi()));
    assert.ok(Number.isFinite(molecule.labuteAsa(false, false)));

    const named = molecule.withName("ethanol").withProperty("source", "binding-test");
    assert.equal(named.nameOrEmpty(), "ethanol");
    assert.equal(named.propertyOrEmpty("source"), "binding-test");
    assert.ok(named.propertyKeys().includes("source"));

    const withSdfField = named.withSdfDataField("ID", "ethanol-1");
    assert.deepEqual([...withSdfField.sdfDataFieldNames()], ["ID"]);
    assert.equal(withSdfField.sdfDataFieldOrEmpty("ID"), "ethanol-1");
    assert.equal(withSdfField.sourceCoordinateDimensionOrEmpty(), "");

    assert.deepEqual([...molecule.atomicNumbers()], [6, 6, 8]);

    const bits = molecule.patternFingerprint();
    assert.ok(bits instanceof Uint32Array);
    assert.ok(bits.length > 0);
    assert.ok([...bits].every((bit) => bit < 2048));
    for (const fingerprint of [
        molecule.morganFingerprint(),
        molecule.atomPairFingerprint(),
        molecule.layeredFingerprint(),
        molecule.topologicalFingerprint(),
        molecule.maccsFingerprint(),
        molecule.avalonFingerprint(),
        molecule.topologicalTorsionFingerprint(),
    ]) {
        assert.ok(fingerprint instanceof Uint32Array);
    }

    const svg = molecule.toSvg(120, 80);
    assert.match(svg, /<svg /);
    assert.match(svg, /width='120px'/);
    assert.match(svg, /height='80px'/);

    const withHydrogens = molecule.withHydrogens();
    assert.ok(withHydrogens.numAtoms() > molecule.numAtoms());
    assert.equal(withHydrogens.withoutHydrogens().numAtoms(), molecule.numAtoms());
    assert.equal(molecule.with2dCoordinates().coordinates2d().length, 6);
    assert.equal(molecule.with2dCoordinates().sourceCoordinateDimensionOrEmpty(), "2D");
    assert.equal(molecule.with3dConformer().numConformers3d(), 1);
    assert.equal(molecule.stereoisomerCount(), "1");
    assert.equal(molecule.hasSubstructMatch("[#8]"), true);

    const binary = molecule.toBinary();
    assert.ok(binary instanceof Uint8Array);
    assert.equal(binding.Molecule.fromBinary(binary).toSmiles(), "CCO");
    assert.equal(molecule.largestFragment().toSmiles(), "CCO");
    assert.equal(molecule.murckoScaffold().numAtoms(), 0);
    assert.equal(molecule.netScaffold().numAtoms(), 0);

    assert.throws(() => binding.Molecule.fromSmiles("["));
    molecule.free();
});
