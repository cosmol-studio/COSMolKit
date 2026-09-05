import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCN.CCO");
console.log("contains oxygen:", molecule.hasSubstructMatch("[#8]"));

const fragments = molecule.fragments();
console.log("fragment SMILES:", fragments.map((fragment) => fragment.toSmiles()));
for (const fragment of fragments) fragment.free();

const largest = molecule.largestFragment();
console.log("largest fragment:", largest.toSmiles());

largest.free();
molecule.free();
