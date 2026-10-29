import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("c1ccccc1O");
console.log(molecule.toSmiles());
console.log(molecule.numAtoms(), molecule.numBonds());

const svg = molecule.toSvg(320, 240);
document.querySelector("#molecule").innerHTML = svg;

molecule.free();
