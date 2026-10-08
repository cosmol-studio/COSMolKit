import {Molecule,PropertyStringError,PropertyValueKind} from 'cosmolkit-generated';const m=Molecule.fromSmiles('CC'),atom:string|null=m.atomPropertyString(0,'x'),bond:string|null=m.bondPropertyString(0,'x');declare const e:PropertyStringError;const kind:PropertyValueKind=e.kind(),domain:string=e.domain,message:string=e.message,cause:Error|null=e.cause;
// @ts-expect-error wasm32 index is a checked number
m.atomPropertyString(1n,'x');
// @ts-expect-error query can return null
const certain:string=m.bondPropertyString(0,'x');
