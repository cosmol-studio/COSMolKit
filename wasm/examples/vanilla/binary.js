import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const source = Molecule.fromSmiles("C[C@H](N)C(=O)O");
const bytes = source.toBinary();
const restored = Molecule.fromBinary(bytes);

console.log({ bytes: bytes.byteLength, smiles: restored.toSmiles() });

restored.free();
source.free();
