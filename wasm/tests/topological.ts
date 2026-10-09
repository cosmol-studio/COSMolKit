import {Molecule,QueryGraph,Fingerprint,TopologicalFingerprintParams,TopologicalFingerprintOutputRequest,TopologicalFingerprintResult,TopologicalFingerprintOutput,TopologicalFingerprintError,fingerprintTopologicalQueryWithParams,fingerprintTopologicalQueryWithOutputWithParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),q=QueryGraph.fromSmarts('[#8]'),p=new TopologicalFingerprintParams(),req=new TopologicalFingerprintOutputRequest(true,true);
const a:Fingerprint=m.fingerprintTopological(),b:Fingerprint=m.fingerprintTopologicalWithParams(p),c:TopologicalFingerprintResult=m.fingerprintTopologicalWithOutput(),d:TopologicalFingerprintResult=m.fingerprintTopologicalWithOutputWithParams(p,req),e:Fingerprint=fingerprintTopologicalQueryWithParams(q,p),f:TopologicalFingerprintResult=fingerprintTopologicalQueryWithOutputWithParams(q,p,req);
const rows:number[][]=f.atomBits(),map:Map<number,number[][]>=f.bitInfo();declare const out:TopologicalFingerprintOutput;const nullable:number[][]|null=out.atomBits,nullableMap:Map<number,number[][]>|null=out.bitInfo;declare const error:TopologicalFingerprintError;const field:string|null=error.field,reason:string|null=error.reason;
// @ts-expect-error typed request is required
m.fingerprintTopologicalWithOutputWithParams(p);
// @ts-expect-error configuration is immutable
p.minPath=2;
// @ts-expect-error queries are not molecules
fingerprintTopologicalQueryWithParams(m,p);
