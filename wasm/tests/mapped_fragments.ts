import {FragmentExtractionParams, type MappedFragment, type Molecule} from 'cosmolkit';

declare const source: Molecule;
declare const names: readonly string[];

const params = new FragmentExtractionParams({sanitize: false});
const parts: readonly MappedFragment[] = source.fragmentsWithAtomMappingWithParams(params);
params.free();

for (const part of [...parts, ...source.fragmentsWithAtomMapping({sanitize: true})]) {
    const molecule: Molecule = part.molecule;
    const indices = part.sourceAtomIndices;
    const fragmentNames: string[] = indices.map(i => names[i]);
    // @ts-expect-error the mapping is read-only
    indices[0] = 0;
    molecule.free();
    part.free();
}
