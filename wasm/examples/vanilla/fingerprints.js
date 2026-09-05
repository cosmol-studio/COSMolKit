import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCOc1ccc2nc(S(N)(=O)=O)sc2c1");
for (const [name, bits] of Object.entries({
  pattern: molecule.patternFingerprint(),
  morgan: molecule.morganFingerprint(),
  atomPair: molecule.atomPairFingerprint(),
  layered: molecule.layeredFingerprint(),
  topological: molecule.topologicalFingerprint(),
  maccs: molecule.maccsFingerprint(),
  avalon: molecule.avalonFingerprint(),
  topologicalTorsion: molecule.topologicalTorsionFingerprint(),
})) {
  console.log(`${name}: ${bits.length} set bits`, bits);
}

molecule.free();
