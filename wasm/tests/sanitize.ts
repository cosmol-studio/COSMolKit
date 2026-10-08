import {Molecule,SanitizeParams,SanitizeOperations,SanitizeStage,SanitizeError,ChemistryProblemReport,ChemistryProblem,ChemistryProblemError,ValenceError,BondOrder,KekulizeParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('C'),p=new SanitizeParams(SanitizeOperations.PROPERTIES),v:Molecule=m.sanitize(),v2:Molecule=m.sanitizeWithParams(p),i:void=m.sanitize_(),i2:void=m.sanitizeWithParams_(p),r:ChemistryProblemReport=m.detectChemistryProblems(),r2:ChemistryProblemReport=m.detectChemistryProblemsWithParams(p);
const problems:ChemistryProblem[]=r.problems,stage:SanitizeStage=problems[0].operation,error:Error=problems[0].error,detail:ChemistryProblemError=problems[0].error.detail,cause:Error=detail.cause;
declare const s:SanitizeError;const ss:SanitizeStage|null=s.stage,bits:number|null=s.bits,unknown:number|null=s.unknownBits;
declare const e:ValenceError;const atom:number|null=e.atom,exp:number|null=e.explicitValence,physical:number|null=e.physicalBonds,atomic:number|null=e.atomicNumber,charge:number|null=e.formalCharge,phase:string|null=e.phase,calculated:number|null=e.calculated,reason:string|null=e.reason,ac:number|null=e.atomCount,neighbor:number|null=e.neighborAtom,bond:number|null=e.bond,bc:number|null=e.bondCount,begin:number|null=e.begin,end:number|null=e.end,value:number|null=e.value,field:string|null=e.field,explicit:number|null=e.explicit,implicit:number|null=e.implicit,nh:number|null=e.neighborHydrogens,order:BondOrder|null=e.order;
// @ts-expect-error immutable report
r.problems=[];
// @ts-expect-error exact sanitize parameter
m.sanitizeWithParams(new KekulizeParams());
// @ts-expect-error mutation returns void
const wrong:Molecule=m.sanitize_();
// @ts-expect-error required params
m.detectChemistryProblemsWithParams();
