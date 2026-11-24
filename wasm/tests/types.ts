// Few witnesses of framework protocols, not a public API coverage inventory.
import {featureSpecs, FeatureSpec} from 'cosmolkit';
import {Reaction} from 'cosmolkit';
const sharedParams = new Reaction().matchParams();
sharedParams.copy().maxMatches = 17;
// @ts-expect-error native shared views grant no writable configuration fields
sharedParams.maxMatches = 17;
for (const value of featureSpecs()) { const typed: FeatureSpec = value; }
import {Molecule, QueryGraph, SubstructMatchParams, TautomerProgress, RingSearchParams, NestedParams, Selection, roundtrip, outcome, bytes, integer, raise, OperationError, RingFindingError, TemplateAttachmentOrderError, Element, ElementParseError} from 'cosmolkit';
const mol = new Molecule('C1CC1');
const params = new RingSearchParams({includeDativeBonds: true});
new NestedParams({child: params}).child.includeDativeBonds = true;
const count: number = mol.symmetrizedSssr(params).numRings;
mol.symmetrizedSssr({includeDativeBonds: true});
mol.symmetrizedSssr();
const matcher = new SubstructMatchParams({extraFinalCheck: (m, indices) => m.numAtoms > count && indices.length > 0});
matcher.extraAtomCheck = (query: QueryGraph, queryAtom: import('cosmolkit').QueryAtom,
    target: Molecule, targetAtom: import('cosmolkit').Atom): boolean => target.numAtoms > 0;
const finalCheck: ((target: Molecule, indices: ReadonlyArray<number>) => boolean) | undefined = matcher.extraFinalCheck;
matcher.extraFinalCheck = undefined;
matcher.extraFinalCheck = null;
matcher.extraFinalCheck = finalCheck;
mol.hasSubstructMatch(QueryGraph.fromSmarts('C'), matcher);
mol.enumerateTautomers({callback: (source, progress) => {
    const saved: TautomerProgress = progress.snapshot();
    return saved.len > 0 && source.numAtoms > 0;
}});
const selection = new Selection({indices: [0, 1]});
selection.indices = [2];
const vector: Uint8Array = bytes(new Uint8Array([1]));
const large: bigint = integer(9007199254740993n);
roundtrip({point: [1, 2, 3], indices: [0], mode: 'preserve'});
const result = outcome({kind: 'completed', count: 1});
if (result.kind === 'completed') { const n: number = result.count; }
// @ts-expect-error optional callback does not make its return type optional
matcher.extraFinalCheck = () => undefined;
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
// @ts-expect-error checked numeric input still has a precise number declaration
Element.fromAtomicNumber('6');
// @ts-expect-error generated setters retain boolean input types
params.includeDativeBonds = 'false';
// @ts-expect-error generated constructor options retain native field types
new RingSearchParams({includeDativeBonds: 'true'});
try { raise(); } catch (error: unknown) {
    if (error instanceof OperationError && error.cause instanceof RingFindingError && error.cause.cause instanceof TemplateAttachmentOrderError) {
        const target: number | undefined = error.cause.cause.details.target;
        const source: Error = error.cause.cause;
    }
}
try { Element.fromSymbol('not-an-element'); } catch (error: unknown) {
    if (error instanceof ElementParseError) {
        const nativeError: Error = error;
        const name: string = error.name;
        const message: string = error.message;
        const input: string = error.details.input;
    }
}
