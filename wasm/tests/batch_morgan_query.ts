import {MoleculeBatch,MorganParams,MorganCallParams,MorganAtomInvariantsGenerator,QueryGraph,BatchQueryParams,Fingerprint} from 'cosmolkit-generated';
import {TopologicalTorsionFingerprintParams} from 'cosmolkit-generated';
const torsions:(Fingerprint|null)[]=MoleculeBatch.fromSmilesList(['CCO']).fingerprintTopologicalTorsionList({},{});
MoleculeBatch.fromSmilesList(['CCO']).fingerprintTopologicalTorsionListWithParams(new TopologicalTorsionFingerprintParams(),new BatchQueryParams());
const atom = MorganAtomInvariantsGenerator.features([QueryGraph.fromSmarts('[#8]')]);
const output: (Fingerprint | null)[] = MoleculeBatch.fromSmilesList(['CCO']).fingerprintMorganListWithGeneratorParams(new MorganParams(2), atom, null, new MorganCallParams(), new BatchQueryParams());
