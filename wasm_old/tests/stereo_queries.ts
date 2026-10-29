import {Molecule, StereoReadError} from 'cosmolkit-generated';
const molecule = Molecule.fromSmiles('CC(O)Cl');
const defaults: Array<[number, string]> = molecule.findChiralCenters();
const explicit: Array<[number, string]> = molecule.findChiralCenters(true);
for (const [atom, label] of explicit) { const index: number = atom, cip: string = label; }
declare const error: StereoReadError;
const kind: string = error.kind, domain: string = error.domain, message: string = error.message, cause: Error | null = error.cause;
// @ts-expect-error selectors require booleans
molecule.findChiralCenters('true');
// @ts-expect-error null does not select the default
molecule.findChiralCenters(null);
