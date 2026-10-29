import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CC(O)C(=O)O");
console.log("stereoisomer count:", molecule.stereoisomerCount());

const threeDimensional = molecule.with3dConformer();
console.log("conformers:", threeDimensional.numConformers3d());
console.log("coordinates:", threeDimensional.coordinates3d(0));

threeDimensional.free();
molecule.free();
