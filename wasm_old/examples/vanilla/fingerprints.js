import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCOc1ccc2nc(S(N)(=O)=O)sc2c1");
try {
  for (const {name, create} of [
    {name: "pattern", create: () => molecule.fingerprintPattern()},
    {name: "morgan", create: () => molecule.fingerprintMorgan()},
    {name: "atomPair", create: () => molecule.fingerprintAtomPair()},
    {name: "layered", create: () => molecule.fingerprintLayered()},
    {name: "topological", create: () => molecule.fingerprintTopological()},
    {name: "maccs", create: () => molecule.fingerprintMaccs()},
    {name: "avalon", create: () => molecule.fingerprintAvalon()},
    {name: "topologicalTorsion", create: () => molecule.fingerprintTopologicalTorsion()},
  ]) {
    const fingerprint = create();
    try {
      const bits = fingerprint.onBits;
      console.log(`${name}: ${bits.length} set bits`, bits);
    } finally {
      fingerprint.free();
    }
  }
} finally {
  molecule.free();
}
