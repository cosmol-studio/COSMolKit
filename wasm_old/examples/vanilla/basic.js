import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("c1ccccc1O");
try {
  console.log(molecule.toSmiles());
  console.log(molecule.numAtoms, molecule.numBonds);
  console.log(molecule.toSvg(320, 240));
} finally {
  molecule.free();
}
