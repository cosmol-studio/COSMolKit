import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { pathToFileURL } from "node:url";

const ETHANOL_SDF = `ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
`;

test("WASM binding preserves runtime values and errors", async () => {
    const modulePath = process.env.COSMOLKIT_WASM_MODULE;
    const wasmPath = process.env.COSMOLKIT_WASM_BINARY;
    assert.ok(modulePath, "COSMOLKIT_WASM_MODULE must point to generated wasm-bindgen JS");
    assert.ok(wasmPath, "COSMOLKIT_WASM_BINARY must point to generated wasm-bindgen binary");

    const binding = await import(pathToFileURL(modulePath).href);
    binding.initSync({ module: readFileSync(wasmPath) });

    const molecule = binding.Molecule.fromSmiles("CCO");
    assert.equal(molecule.numAtoms(), 3);
    assert.equal(molecule.numBonds(), 2);
    assert.deepEqual([...molecule.atomicNumbers()], [6, 6, 8]);
    assert.equal(binding.Molecule.fromSmilesWithSanitize("CCO", false).numAtoms(), 3);
    assert.equal(binding.Molecule.fromSdf(ETHANOL_SDF).numAtoms(), 3);

    assert.equal(molecule.molecularWeight(), 46.069);
    assert.ok(Math.abs(molecule.exactMolecularWeight() - 46.041864812) < 1e-9);
    assert.equal(molecule.molecularFormula(), "C2H6O");

    const named = molecule.withName("ethanol").withProperty("source", "binding-test");
    assert.equal(named.nameOrEmpty(), "ethanol");
    assert.equal(named.propertyOrEmpty("source"), "binding-test");
    assert.equal(named.propertyOrEmpty("absent"), "");
    assert.ok([...named.propertyKeys()].includes("source"));

    const withSdfField = named.withSdfDataField("ID", "ethanol-1");
    assert.deepEqual([...withSdfField.sdfDataFieldNames()], ["ID"]);
    assert.equal(withSdfField.sdfDataFieldOrEmpty("ID"), "ethanol-1");
    assert.equal(withSdfField.sdfDataFieldOrEmpty("absent"), "");

    const withHydrogens = molecule.withHydrogens();
    assert.equal(withHydrogens.numAtoms(), 9);
    assert.equal(withHydrogens.withoutHydrogens().numAtoms(), 3);

    assert.equal(molecule.numConformers3d(), 0);
    assert.equal(molecule.coordinates2d().length, 0);
    const drawn = molecule.with2dCoordinates();
    assert.equal(drawn.coordinates2d().length, 6);
    assert.equal(drawn.numAtoms(), 3);

    const benzene = binding.Molecule.fromSmiles("c1ccccc1");
    benzene.withKekulizedBonds();
    benzene.withAssignedAromaticity();
    molecule.sanitize();
    molecule.withAssignedValence();
    molecule.withAssignedRings();
    molecule.withAssignedRingFamilies();
    molecule.withAssignedRadicals();

    assert.throws(() => binding.Molecule.fromSmiles("["));
    assert.throws(() => binding.Molecule.fromSdf("not an SDF record"));
    molecule.free();
});
