import {Molecule,MoleculeHashError,CipRankError} from 'cosmolkit-generated';const m=Molecule.fromSmiles('CCO'),a:bigint=m.molecularHash(),b:bigint=m.molecularHashWithRanks([0,1,2]),c:bigint=m.molecularHashWithRanks(new Uint32Array([0,1,2]));declare const e:MoleculeHashError,p:CipRankError;const kind:string=e.kind,cause:Error|null=e.cause,actual:number|null=e.actual,atoms:number|null=e.atomCount,map:number|null=p.mapNumber,value:number|bigint|null=p.value,field:string|null=p.field;
const scaffold: Molecule = m.murckoScaffold(), net: Molecule = m.netScaffold(), decomposition: Molecule = m.murckoDecompose();
// @ts-expect-error unsigned ranks are numbers
m.molecularHashWithRanks([0n,1n,2n]);
// @ts-expect-error u64 is bigint
const lossy:number=m.molecularHash();
