import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCOc1ccc2nc(S(N)(=O)=O)sc2c1");
for (const [name, bits] of Object.entries({
  pattern: molecule.fingerprintPattern(),
  morgan: molecule.fingerprintMorgan(),
  atomPair: molecule.fingerprintAtomPair(),
  layered: molecule.fingerprintLayered(),
  topological: molecule.fingerprintTopological(),
  maccs: molecule.fingerprintMaccs(),
  avalon: molecule.avalonFingerprint(),
  topologicalTorsion: molecule.fingerprintTopologicalTorsion(),
})) {
  console.log(`${name}: ${bits.length} set bits`, bits);
}

molecule.free();
