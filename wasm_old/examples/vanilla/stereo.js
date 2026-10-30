import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CC(O)C(=O)O");
try {
  console.log("stereoisomer count:", molecule.stereoisomerCount());
  const hydrogenated = molecule.withHydrogens();
  try {
    const threeDimensional = hydrogenated.with3dConformer({randomSeed: 42});
    try {
      console.log("conformers:", threeDimensional.num3dConformers);
      console.log("coordinates:", threeDimensional.coordinates3d(0));
    } finally {
      threeDimensional.free();
    }
  } finally {
    hydrogenated.free();
  }
} finally {
  molecule.free();
}
