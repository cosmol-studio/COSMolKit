import {MoleculeBatch,MorganParams,MorganInvariants,MorganFingerprintParams,MorganCallParams,MorganAtomInvariantsGenerator,MorganBondInvariantsGenerator,MorganReadError,QueryGraph,BatchQueryParams,Fingerprint,BatchFingerprintOutput} from 'cosmolkit-generated';
const generator=new MorganParams(2),options=new MorganFingerprintParams(generator,[0],[],null,null,-1,MorganInvariants.features()),call=MorganCallParams.new(null,[],null,null,-1),query=new BatchQueryParams(),atom=MorganAtomInvariantsGenerator.features([QueryGraph.fromSmarts('[#8]')]),bond=MorganBondInvariantsGenerator.new(true,false),batch=MoleculeBatch.fromSmilesList(['CCO']);
const bits:(Fingerprint|null)[][]=[batch.fingerprintMorganList(),batch.fingerprintMorganListWithParams(options,query),batch.fingerprintMorganListWithGeneratorParams(generator,atom,bond,call,query)];
const outputs:(BatchFingerprintOutput|null)[][]=[batch.fingerprintMorganWithOutputList(),batch.fingerprintMorganWithOutputListWithParams(options,true,query),batch.fingerprintMorganWithOutputListWithGeneratorParams(generator,null,null,new MorganCallParams(),false,query)];
const p:MorganParams=options.generator,r:number=p.radius,roots:number[]|null=call.fromAtoms,flags:boolean=bond.useBondTypes();declare const error:MorganReadError;const kind:string=error.kind;
// @ts-expect-error readonly generator
p.radius=9;
// @ts-expect-error query graph values required
MorganAtomInvariantsGenerator.features(['O']);
// @ts-expect-error typed boolean collection policy
batch.fingerprintMorganWithOutputListWithParams(options,'false',query);
void [bits,outputs,r,roots,flags,kind];
