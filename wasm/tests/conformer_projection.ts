import {Molecule,EmbedParams,EmbedMoleculeResult,EmbedMultipleConfsResult,ConformerError,ConformerRunError} from 'cosmolkit-generated';
const p=new EmbedParams(0,1,42),factories:EmbedParams[]=[EmbedParams.new(),EmbedParams.dg(),EmbedParams.kdg(),EmbedParams.etdg(),EmbedParams.etdgV2(),EmbedParams.etkdg(),EmbedParams.etkdgV2(),EmbedParams.etkdgV3(),EmbedParams.srEtkdgV3()],m=Molecule.fromSmiles('C');const cm:Map<number,number[]>|null=p.coordMap(),cpci:Map<readonly [number,number],number>|null=p.cpci(),failures:Uint32Array=p.failures(),json:string=p.toJson(),changed:EmbedParams=p.withJson(json),bounds:number[][]=m.dgBoundsMatrix();
const r0:Molecule=m.with3dConformer();
const r1:Molecule=m.with3dConformerWithParams(p);
const r2:void=m.embed3dConformer();
const r3:void=m.embed3dConformerWithParams(p);
const r4:EmbedMoleculeResult=m.with3dConformerResult();
const r5:EmbedMoleculeResult=m.with3dConformerResultWithParams(p);
const r6:EmbedMoleculeResult=m.embed3dConformerResult();
const r7:EmbedMoleculeResult=m.embed3dConformerResultWithParams(p);
const r8:Molecule=m.with3dConformers(2);
const r9:Molecule=m.with3dConformersWithParams(2,p);
const r10:void=m.embed3dConformers(2);
const r11:void=m.embed3dConformersWithParams(2,p);
const r12:EmbedMultipleConfsResult=m.with3dConformersResult(2);
const r13:EmbedMultipleConfsResult=m.with3dConformersResultWithParams(2,p);
const r14:EmbedMultipleConfsResult=m.embed3dConformersResult(2);
const r15:EmbedMultipleConfsResult=m.embed3dConformersResultWithParams(2,p);
declare const single:EmbedMoleculeResult;declare const multiple:EmbedMultipleConfsResult;const returned:Molecule=single.molecule(),options:EmbedParams=single.params(),ok:boolean=single.ok(),id:number=single.confId(),ids:Int32Array=multiple.confIds(),requested:number=multiple.requestedNumConfs(),generated:number=multiple.generatedCount();declare const err:ConformerError;declare const run:ConformerRunError;const detail:string|null=err.detail,cause:Error|null=run.cause;
// @ts-expect-error checked number argument
m.with3dConformersWithParams('2',p);
// @ts-expect-error readonly parameter methods
p.randomSeed=42;
void [factories,cm,cpci,failures,changed,bounds,returned,options,ok,id,ids,requested,generated,detail,cause];
