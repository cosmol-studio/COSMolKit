import { Molecule, MoleculeBatch, MolBlockWriteParams, SdfRecord } from "cosmolkit-generated";

const molecule: Molecule = Molecule.fromSmiles("CCO");
const params = new MolBlockWriteParams(undefined, undefined, false);
const text: string = molecule.toSdfWithParams(params);
const record: SdfRecord = SdfRecord.fromSdf(text);
const restored: Molecule = record.molecule();
const batch: MoleculeBatch = MoleculeBatch.fromSmilesList(["CCO"]);
batch.toSmilesList();

// Native archives are absent even from the full WASM distribution.
// @ts-expect-error platform-excluded native serialization
molecule.toBinary();
// @ts-expect-error platform-excluded native serialization
Molecule.fromBinary(new Uint8Array());

restored.free(); record.free(); params.free(); molecule.free(); batch.free();
