import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
try {
  const crippen = molecule.crippenDescriptors();
  try {
    console.table({
      formula: molecule.molecularFormula(),
      molecularWeight: molecule.molecularWeight(),
      exactMolecularWeight: molecule.exactMolecularWeight(),
      logP: crippen.logp,
      tpsa: molecule.tpsa(),
      rings: molecule.numRings(),
      heavyAtoms: molecule.numHeavyAtoms(),
    });
  } finally {
    crippen.free();
  }
  const withCoordinates = molecule.with2dCoordinates();
  try {
    console.log("2D coordinate values:", withCoordinates.coordinates2d);
  } finally {
    withCoordinates.free();
  }
} finally {
  molecule.free();
}
