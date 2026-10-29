import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
console.table({
  formula: molecule.formula(),
  molecularWeight: molecule.molecularWeight(),
  exactMolecularWeight: molecule.exactMolecularWeight(),
  logP: molecule.crippenLogP(),
  tpsa: molecule.tpsa(),
  rings: molecule.numRings(),
  heavyAtoms: molecule.numHeavyAtoms(),
});

const withCoordinates = molecule.with2dCoordinates();
console.log("2D coordinate values:", withCoordinates.coordinates2d());

withCoordinates.free();
molecule.free();
