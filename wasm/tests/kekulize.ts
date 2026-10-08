import {Molecule,KekulizeParams,KekulizeError,RemoveHsParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('c1ccccc1'),p=new KekulizeParams(false,true,4294967295);
const a:Molecule=m.withKekulizedBonds(),b:Molecule=m.withKekulizedBondsWithParams(p);
const c:void=m.kekulizeBonds(),d:void=m.kekulizeBondsWithParams(p);
const mark:boolean=p.markAtomsBonds,canonical:boolean=p.canonical,max:number=p.maxBacktracks;
declare const e:KekulizeError;
const domain:string=e.domain,kind:string=e.kind,message:string=e.message,cause:Error|null=e.cause;
const expected:number|null=e.expected,actual:number|null=e.actual,atom:number|null=e.atom,atomCount:number|null=e.atomCount;
const field:string|null=e.field,begin:number|null=e.begin,end:number|null=e.end,questions:number|null=e.questions,bitWidth:number|null=e.bitWidth;
const problemAtoms:number[]|null=e.problemAtoms,before:number|null=e.before,after:number|null=e.after,bond:number|null=e.bond,reason:string|null=e.reason;
// @ts-expect-error exact parameter class
m.withKekulizedBondsWithParams(new RemoveHsParams());
// @ts-expect-error in-place returns void
const result:Molecule=m.kekulizeBonds();
// @ts-expect-error u32 represented by number
new KekulizeParams(true,true,1n);
// @ts-expect-error parameters required
m.kekulizeBondsWithParams();
