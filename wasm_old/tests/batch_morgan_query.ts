import {MoleculeBatch,MorganGeneratorParams,MorganCallParams,MorganAtomInvariantsGenerator,QueryGraph,BatchQueryParams,Fingerprint} from 'cosmolkit-generated';
import {TopologicalTorsionFingerprintParams} from 'cosmolkit-generated';
const torsions:(Fingerprint|null)[]=MoleculeBatch.fromSmilesList(['CCO']).fingerprintTopologicalTorsion({},{});
MoleculeBatch.fromSmilesList(['CCO']).fingerprintTopologicalTorsionWithParams(new TopologicalTorsionFingerprintParams(),new BatchQueryParams());
const atom = MorganAtomInvariantsGenerator.features([QueryGraph.fromSmarts('[#8]')]);
const output: (Fingerprint | null)[] = MoleculeBatch.fromSmilesList(['CCO']).fingerprintMorganWithGeneratorParams(new MorganGeneratorParams(2), atom, null, new MorganCallParams(), new BatchQueryParams());
