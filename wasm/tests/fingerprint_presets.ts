import { Molecule, Fingerprint, MorganFingerprintGenerator, MorganAtomInvariantsGenerator, MorganParams, TopologicalFingerprintParams, TopologicalFingerprintOutputRequest, LayeredFingerprintParams, PatternFingerprintParams } from "cosmolkit-generated";
const molecule = Molecule.fromSmiles("CCO");
const generator = new MorganFingerprintGenerator(new MorganParams(0), MorganAtomInvariantsGenerator.features([]));
const morgan: Fingerprint = molecule.morganFingerprintWithGenerator(generator);
const topological = molecule.topologicalFingerprintWithOutputWithParams(new TopologicalFingerprintParams(), new TopologicalFingerprintOutputRequest(true, true));
const paths: Map<number, number[][]> = topological.bitInfo();
const layered: Fingerprint = molecule.layeredFingerprintWithParams(new LayeredFingerprintParams());
const pattern: Fingerprint = molecule.patternFingerprintWithParams(new PatternFingerprintParams());
molecule.free(); generator.free(); morgan.free(); layered.free(); pattern.free(); topological.free();
