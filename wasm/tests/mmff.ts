import {Molecule,MmffEvaluationParams,MmffOptimizationParams,MmffConformerOptimizationParams,MmffEnergyGradient,MmffOptimizeMoleculeResult,MmffOptimizeMoleculeConfsResult,MmffOptimizeMoleculeConfResult,MmffOptimizationError} from 'cosmolkit-generated';const m=Molecule.fromSdf(''),ep=new MmffEvaluationParams(),sp=new MmffOptimizationParams('MMFF94s',2,100,-1,false),cp=new MmffConformerOptimizationParams();const a:MmffEnergyGradient|null=m.mmffEnergyGradient(),b:MmffEnergyGradient|null=m.mmffEnergyGradientWithParams(ep),s:MmffOptimizeMoleculeResult=m.withMmffOptimized(),t:MmffOptimizeMoleculeResult=m.withMmffOptimizedWithParams(sp),c:MmffOptimizeMoleculeConfsResult=m.withMmffOptimizedConfs(),d:MmffOptimizeMoleculeConfsResult=m.withMmffOptimizedConfsWithParams(cp);if(a){const energy:number=a.energy(),gradient:number[]=a.gradient();}const out:Molecule=s.molecule(),status:number=t.statusCode(),more:boolean=t.needsMore(),rows:MmffOptimizeMoleculeConfResult[]=c.conformerResults(),other:Molecule=d.molecule(),re:number=rows[0].energy(),rm:boolean=rows[0].needsMore(),rs:number=rows[0].statusCode(),selector:number|null=ep.conformerId;declare const e:MmffOptimizationError;const cause:Error|null=e.cause;
// @ts-expect-error nullable evaluation
a.energy();
// @ts-expect-error typed parameter
m.mmffEnergyGradientWithParams('bad');
// @ts-expect-error immutable
sp.maxIterations=1;
// @ts-expect-error strict scalar
new MmffEvaluationParams(94);
