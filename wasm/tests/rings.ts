import {Molecule,RingSearchParams,KekulizeParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('C1CCCCC1'),p=new RingSearchParams(true,true),d:boolean=p.includeDativeBonds,h:boolean=p.includeHydrogenBonds;
const a:Molecule=m.withAssignedRings(),b:void=m.assignRings(),c:Molecule=m.withAssignedRingFamilies(),e:Molecule=m.withAssignedRingFamiliesWithParams(p),f:void=m.assignRingFamilies(),g:void=m.assignRingFamiliesWithParams(p);
// @ts-expect-error exact parameter class
m.assignRingFamiliesWithParams(new KekulizeParams());
// @ts-expect-error parameters required
m.withAssignedRingFamiliesWithParams();
// @ts-expect-error immutable parameters
p.includeDativeBonds=false;
// @ts-expect-error mutation returns void
const out:Molecule=m.assignRings();
