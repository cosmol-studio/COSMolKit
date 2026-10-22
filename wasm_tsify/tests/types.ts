// Few witnesses of framework protocols, not a public API coverage inventory.
import {Molecule, RingSearchParams, NestedParams, Selection, roundtrip, outcome, bytes, integer, failure_details} from '../pkg/cosmolkit_tsify_probe';
const mol = new Molecule('C1CC1');
const params = new RingSearchParams();
new NestedParams(params).child.include_dative_bonds = true;
const count: number = mol.num_rings(params);
mol.inspect(m => m.num_atoms() > count);
const selection = new Selection([0, 1]);
selection.indices = [2];
const vector: Uint8Array = bytes(new Uint8Array([1]));
const large: bigint = integer(9007199254740993n);
roundtrip({point: [1, 2, 3], indices: [0], mode: 'preserve'});
const result = outcome({kind: 'completed', count: 1});
if (result.kind === 'completed') { const n: number = result.count; }
// @ts-expect-error callback return type is not arbitrary
mol.inspect(() => 3);
// @ts-expect-error frozen container getter cannot be mutated
selection.indices.push(1);
// @ts-expect-error coordinate dimension is fixed
roundtrip({point: [1, 2], indices: [], mode: 'preserve'});
// @ts-expect-error enum spelling is exact
roundtrip({point: [1, 2, 3], indices: [], mode: 'old'});
// @ts-expect-error payload enum is discriminated
outcome({kind: 'completed', message: 'not a count'});
// @ts-expect-error 64-bit integers are not JS numbers
integer(1);
const error = failure_details(new Error());
if (error.kind === 'Syntax') { const offset: number = error.offset; }
