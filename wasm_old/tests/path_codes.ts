import {Molecule,AtomPairsParameters,AtomCodeExplanation,AtomCodeExplanationError,TopologicalTorsionPathScoreError,AtomPairAtomCodeResult,explainPathScore} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCCC');const s:bigint=m.topologicalTorsionPathScore([0,1,2,3],4);m.topologicalTorsionPathScore(new Uint32Array([0]),1,null);m.topologicalTorsionPathScore([0],1,new Uint32Array([33,34,34,33]));const d:[string,number,number][]=explainPathScore(s);explainPathScore(s,0);
const p=AtomPairsParameters;const v:string=p.version();const types:Uint32Array=p.atomTypes();const nums:number[]=[p.numTypeBits(),p.numPiBits(),p.numBranchBits(),p.numChiralBits(),p.codeSize(),p.numPathBits(),p.maxPathLength(),p.numAtomPairFingerprintBits()];const e=AtomCodeExplanation.fromCode(33n,-1n,true);const symbol:string=e.symbol();const branch:number=e.branchCount();const pi:number=e.piElectrons();const chirality:string|null=e.chirality();const r:AtomPairAtomCodeResult=m.withAtomPairAtomCode(0);m.withAtomPairAtomCode(0,0,true,false);const n:Molecule=r.molecule;const code:number=r.code;
declare const error:TopologicalTorsionPathScoreError;const actual:number|null=error.actual;const required:number|null=error.required;const index:number|null=error.index;const atoms:number|null=error.atomCount;const errcode:number|null=error.code;const sub:number|null=error.subtract;const cause:Error|null=error.cause;declare const err:AtomCodeExplanationError;const small:number=err.code;
// @ts-expect-error no floating integer replacement for exact bigint
explainPathScore(1);
// @ts-expect-error signed 64-bit transport is bigint
AtomCodeExplanation.fromCode(33n,0);
// @ts-expect-error result is immutable
r.code=2;
// @ts-expect-error size required
m.topologicalTorsionPathScore([0]);
