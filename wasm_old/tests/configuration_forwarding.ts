import {Molecule, MoleculeBatch, SmilesParseParams, SmilesWriteParams, BatchParams, BatchQueryParams, MorganFingerprintParams, MorganGeneratorParams, AtomPairFingerprintParams, TopologicalTorsionFingerprintParams, Reaction, ReactionWriteParams, AlignmentParameters, MolBlockWriteParams} from 'cosmolkit-generated';
const molecule=Molecule.fromSmiles('CCO');
molecule.fingerprintMorgan(new MorganFingerprintParams());
molecule.fingerprintMorgan({generator:new MorganGeneratorParams()});
molecule.fingerprintAtomPair(new AtomPairFingerprintParams());
molecule.fingerprintTopologicalTorsion(new TopologicalTorsionFingerprintParams());
molecule.toMol(new MolBlockWriteParams());
molecule.toSdf({includeStereo:true});
molecule.withUffOptimized({maxIterations:2});
molecule.withMmffOptimized({maxIterations:2});
molecule.alignmentTransformTo(molecule,new AlignmentParameters());
molecule.alignmentTransformTo(molecule,{reflect:false});
const batch=MoleculeBatch.fromSmilesList(['CCO'],new SmilesParseParams(),new BatchParams());
batch.toSmiles(new SmilesWriteParams(),new BatchQueryParams());
batch.toSmiles({allBondsExplicit:true},{progressCallback:()=>{}});
batch.withHydrogens({},{});
batch.fingerprintMorgan({},{});
Reaction.fromSmirks('[C:1]>>[C:1]',{}).toSmirks(new ReactionWriteParams());
// @ts-expect-error checked optional configuration
molecule.withUffOptimized({maxIterations:'2'});
// @ts-expect-error checked batch parsing configuration
MoleculeBatch.fromSmilesList(['CCO'],{sanitize:'true'});
