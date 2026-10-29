import {MoleculeBatch,MorganGeneratorParams,MorganInvariants,MorganFingerprintParams,MorganCallParams,MorganAtomInvariantsGenerator,MorganBondInvariantsGenerator,MorganReadError,QueryGraph,BatchQueryParams,Fingerprint,BatchFingerprintOutput} from 'cosmolkit-generated';
const generator=new MorganGeneratorParams(2),options=new MorganFingerprintParams(generator,[0],[],null,null,MorganInvariants.features()),call=MorganCallParams.new(null,[],null,null),query=new BatchQueryParams(),atom=MorganAtomInvariantsGenerator.features([]),bond=MorganBondInvariantsGenerator.new(true,false),batch=MoleculeBatch.fromSmilesList(['CCO']);
const bits:(Fingerprint|null)[][]=[batch.fingerprintMorgan(),batch.fingerprintMorganWithParams(options,query),batch.fingerprintMorganWithGeneratorParams(generator,atom,bond,call,query)];
const outputs:(BatchFingerprintOutput|null)[][]=[batch.fingerprintMorganWithOutput(),batch.fingerprintMorganWithOutputWithParams(options,true,query),batch.fingerprintMorganWithOutputWithGeneratorParams(generator,null,null,new MorganCallParams(),false,query)];
const p:MorganGeneratorParams=options.generator,r:number=p.radius,roots:number[]|null=call.fromAtoms,flags:boolean=bond.useBondTypes;declare const error:MorganReadError;const kind:string=error.kind;
p.radius=9;
// @ts-expect-error query graph values required
MorganAtomInvariantsGenerator.features(['O']);
// @ts-expect-error typed boolean collection policy
batch.fingerprintMorganWithOutputWithParams(options,'false',query);
void [bits,outputs,r,roots,flags,kind];
