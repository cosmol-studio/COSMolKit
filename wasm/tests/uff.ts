import {Molecule,UffEvaluationParams,UffOptimizationParams,UffConformerOptimizationParams,UffEnergyGradient,UffOptimizationResult,UffConformerOptimizationResult,UffConformerResult,UffOptimizationError,UffOptimizationErrorKind} from 'cosmolkit-generated';
const m=Molecule.fromSdf(''),ep=new UffEvaluationParams(),sp=new UffOptimizationParams(2,10,-1,false),cp=new UffConformerOptimizationParams();const a:UffEnergyGradient=m.uffEnergyGradient(),b:UffEnergyGradient=m.uffEnergyGradientWithParams(ep),s:UffOptimizationResult=m.withUffOptimized(),t:UffOptimizationResult=m.withUffOptimizedWithParams(sp),c:UffConformerOptimizationResult=m.withUffOptimizedConfs(),d:UffConformerOptimizationResult=m.withUffOptimizedConfsWithParams(cp);const energy:number=a.energy(),gradient:number[]=b.gradient(),out:Molecule=s.molecule(),status:number=t.statusCode(),more:boolean=t.needsMore(),se:number=s.energy(),rows:UffConformerResult[]=c.conformerResults(),other:Molecule=d.molecule(),id:number=rows[0].conformerId(),re:number=rows[0].energy(),rm:boolean=rows[0].needsMore(),rs:number=rows[0].statusCode(),selector:number|null=ep.conformerId;declare const e:UffOptimizationError;const k:UffOptimizationErrorKind=e.kind(),v:string=k.variant,requested:number|null=k.requested,same:boolean=k.equals(k),cause:Error|null=e.cause;
// @ts-expect-error typed parameter
m.uffEnergyGradientWithParams('bad');
// @ts-expect-error immutable
sp.maxIterations=1;
// @ts-expect-error strict scalar
new UffEvaluationParams('10');
